function [Coords_F,Connectivity_F,Lattice_vectors]=Generate_Lattice(THETA,L1,L2,L3)

[n1,c1]=lattice_data(THETA.lat(1));
[n2,c2]=lattice_data(THETA.lat(2));
[n3,c3]=lattice_data(THETA.lat(3));
% truss_plotting(n3,c3)

t1=1+THETA.tess(1);
t2=1+THETA.tess(2);
t3=1+THETA.tess(3);

[Rotax1,Rotan1]=ang_ax(THETA.Rot1);
[Rotax2,Rotan2]=ang_ax(THETA.Rot2);
% Rotan2=-0.35;
% Rotan1=0.57;
lattice_vecs = {L1, L2, L3};

% load individual lattices
lattice_1 = tessellate(n1, c1, t1, lattice_vecs);
lattice_1 = remove_redundancy(lattice_1);

lattice_2 = tessellate(n2, c2, t2, lattice_vecs);
lattice_2 = remove_redundancy(lattice_2);

lattice_3 = tessellate(n3, c3, t3, lattice_vecs);
lattice_3 = remove_redundancy(lattice_3);
% truss_plotting(lattice_3.nodes,lattice_3.connectivity)

% add all the lattices together
list = {lattice_1, lattice_2, lattice_3};
lattice = add_lattices(list);
lattice = remove_redundancy(lattice);

% process the combined lattice for intersections
lattice = process_lattice(lattice);
lattice = remove_redundancy(lattice);

% Run remove_redundancy() again (just to be sure)
lattice = remove_redundancy(lattice);

Coords=lattice.nodes;

U=diag(THETA.stretch1);
V=diag(THETA.stretch2);
R1=axang2rotm([Rotax1, Rotan1]);
% R1=rotate_arbt_axis(Rotax1, Rotan1);
R2=axang2rotm([Rotax2, Rotan2]);
% R2=rotate_arbt_axis(Rotax2, Rotan2);

Coords_F=R2*V*R1*U*Coords';
Coords_F=Coords_F';
Lattice_vectors=R2*V*R1*U*[L1;L2;L3]';
Lattice_vectors=Lattice_vectors';
Connectivity_F=lattice.connectivity;
% truss_plotting(Coords_F,Connectivity_F)
end


function [lattice] = tessellate(nodes, connectivity, tess, lattice_vecs)

% nodes, connectivity, tessellation, lattice vectors 1, 2, & 3
% Note: lattice vectors must be for the ** untessellated ** cell 

if tess==1
    % do nothing
else

    vec1 = lattice_vecs{1}/tess;
    vec2 = lattice_vecs{2}/tess;
    vec3 = lattice_vecs{3}/tess;

    n_unit = nodes/tess;
    c_unit = connectivity;

    nodes=[];
    connectivity=[];

    for i=1:tess
        for j=1:tess
            for k=1:tess
                
                %order of next two lines is critical:
                connectivity = [connectivity; c_unit + size(nodes,1)];
                nodes = [nodes; (n_unit + (i-1)*vec1 + (j-1)*vec2 + (k-1)*vec3)];
            end
        end
    end

end

lattice = get_lattice(nodes, connectivity);

end

function [lattice] = get_lattice(nodes, connectivity)
lattice.nodes = nodes;
lattice.connectivity = connectivity;
end

function [new_lattice] = remove_redundancy(lattice)

nodes = lattice.nodes;
connectivity = lattice.connectivity;

cutoff_radius = 0.01;

% identify and map duplicate nodes to node with lowest id
nodes_map = zeros(size(nodes,1),1);

idx = rgsearch(nodes,nodes,cutoff_radius);

for i=1:size(nodes,1)
    nodes_map(i)=min(idx{i});
end

for c=1:size(connectivity,1)
   for k=1:2
       connectivity(c,k) = nodes_map(connectivity(c,k));
   end
end

% % remove connectivity with duplicate nodes (zero length)
idx = (connectivity(:,1) == connectivity(:,2));
connectivity(idx,:) = [];

% remove duplicate connectivity
% first sort connectivity with lower node id on left
for c = 1:size(connectivity,1)
   connectivity(c,:)=sort(connectivity(c,:));
end
[connectivity,~,~]=unique(connectivity,'rows');

% remove unconnected nodes
idx = unique(connectivity);
nodes_map = zeros(1,size(nodes,1));
nodes_map(idx) = 1:size(idx,1);
for c=1:size(connectivity,1)
    for k=1:2
        connectivity(c,k) = nodes_map(connectivity(c,k));
    end
end
nodes(nodes_map==0,:)=[];

new_lattice = get_lattice(nodes, connectivity);

end

function [lattice] = add_lattices(list)

nodes = [];
connectivity = [];

for i = 1:length(list)
    %order of next two lines is critical:
    connectivity = [connectivity; list{i}.connectivity + size(nodes,1)];
    nodes = [nodes; list{i}.nodes];
end

lattice = get_lattice(nodes, connectivity);

end

function [lattice] = process_lattice(lattice)

lattice = add_intersectional_nodes(lattice);
% Do NOT run remove_redundancy() here
lattice = split_connections(lattice);

end

function [new_lattice] = add_intersectional_nodes(lattice)
% Aim is to identify intersection points that are potential new nodes.
% When two elements are parallel/partially coincident/totally coincident,
% one of the existing nodes will automatically be considered as potential
% site for splitting the elements in the next step, and hence nothing needs
% to be done here for those cases. When the two lines are skewed, then
% nothing needs to be done either. However, when the lines are not parallel
% or skewed, two cases arise - (a) the intersection point lies in the
% middle/inside of the both elements, or (b) intersection point coinicides
% with an end point of one or more elements. In the case of (b), the
% node/end-point itself will serve as a candidate for splitting elements
% in the next step, and so nothing needs to be done here. Only case (a)
% needs to be taken care of here.

nodes = lattice.nodes;
connectivity = lattice.connectivity;

nodes_to_add = [];

for j = 1 : size(connectivity,1)
    
    for k = (j+1) : size(connectivity,1)
        
        p1 = nodes(connectivity(j,1),:);
        p2 = nodes(connectivity(j,2),:);

        p3 = nodes(connectivity(k,1),:);
        p4 = nodes(connectivity(k,2),:);

        if check_parallel_lines(p1, p2, p3, p4)
            % lines are parallel or coincident; nothing to do
            continue; 
        elseif check_skew_lines(p1, p2, p3, p4)
            % lines are parallel or coincident; nothing to do
            continue;
        else
            % lines are intersecting
            [px, check_inside_segments] = intersection_point(p1, p2, p3, p4);
            
            if check_inside_segments
                nodes_to_add = [nodes_to_add; px];
            end
        end

    end % for loop: k

end % for loop: j

% if nodes_to_add is NOT empty
if ~isempty(nodes_to_add)
    % find unique nodes in nodes_to_add
    nodes_to_add = find_unique_nodes(nodes_to_add);
    
    % remove ones in nodes_to_add that already exist in nodes
    nodes_to_add = find_new_nodes(nodes_to_add, nodes);
    
    % add intersectional nodes
    nodes = [nodes; nodes_to_add];
end

new_lattice = get_lattice(nodes,connectivity);

end


function [val] = check_parallel_lines(p1, p2, p3, p4)

vecA = p1-p2;
vecA = vecA/norm(vecA);

vecB = p3-p4;
vecB = vecB/norm(vecB);

delta = abs(dot(vecA, vecB)) - 1.0;

val = (abs(delta) < 1e-8);
end


function [val] = check_skew_lines(p1, p2, p3, p4)
% https://mathworld.wolfram.com/SkewLines.html

A = [p1; p2; p3; p4];

A(:,4) = 1;

d = abs(det(A));

val = (abs(d) > 1e-8);
end


function [nodes] = find_unique_nodes(nodes)
% identify and compare duplicate nodes to node with lowest id
cutoff_radius = 0.01;
nodes_flag = (zeros(size(nodes,1),1)==0);
idx = rgsearch(nodes, nodes, cutoff_radius);
for i=1:size(nodes,1)
    nodes_flag(i) = (i == min(idx{i}));
end
nodes = nodes(nodes_flag,:);
end

function [nodes] = find_new_nodes(nodes, reference)
%finds nodes that dont exist in reference
cutoff_radius = 0.01;
idx = rgsearch(reference, nodes, cutoff_radius);
nodes_flag = (zeros(size(nodes,1),1)==0);
for i=1:size(nodes,1)
    nodes_flag(i) = isempty(idx{i});
end
nodes = nodes(nodes_flag,:);
end

function [px, check_inside_segments] = intersection_point(p1, p2, p3, p4)

% Line1: from P1 to P2
[x1, y1, z1] = get_coords_from_point(p1);
[x2, y2, z2] = get_coords_from_point(p2);

% Line2: from P3 to P4
[x3, y3, z3] = get_coords_from_point(p3);
[x4, y4, z4] = get_coords_from_point(p4);

% setup equations (see accompanying mathematica derivations)
A = zeros(2,2);
b = zeros(2,1);

A(1,1) = power(x1,2) - 2*x1*x2 + power(x2,2) + power(y1,2) - 2*y1*y2 + power(y2,2) + power(z1,2) - 2*z1*z2 + power(z2,2);

A(1,2) = -(x1*x3) + x2*x3 + x1*x4 - x2*x4 - y1*y3 + y2*y3 + y1*y4 - y2*y4 - z1*z3 + z2*z3 + z1*z4 - z2*z4;

A(2,1) = x1*x3 - x2*x3 - x1*x4 + x2*x4 + y1*y3 - y2*y3 - y1*y4 + y2*y4 + z1*z3 - z2*z3 - z1*z4 + z2*z4;

A(2,2) = -power(x3,2) + 2*x3*x4 - power(x4,2) - power(y3,2) + 2*y3*y4 - power(y4,2) - power(z3,2) + 2*z3*z4 - power(z4,2);

b(1,1) = -((x1 - x2)*(-x1 + x3)) - (y1 - y2)*(-y1 + y3) - (z1 - z2)*(-z1 + z3);

b(2,1) = -((-x1 + x3)*(x3 - x4)) - (-y1 + y3)*(y3 - y4) - (-z1 + z3)*(z3 - z4);

% check for singularity
if abs(det(A))<1e-8
    error('Equations are singular for intersection point!')
end

% solve equations
sol = A\b;

% check for finite solution
if ~all(isfinite(sol))
    error('Non-finite solution found for intersection point!')
end

% compute intersection points for both lines
mua = sol(1);
mub = sol(2);
pa = p1 + mua * (p2-p1);
pb = p3 + mub * (p4-p3);

% sanity check: pa == pb
delta = norm(pa-pb);
if (delta > 1e-8)
    error('Incorrect solution found for intersection point!')
end

% intersection point: px = pa = pb
px = pa;


check_inside_segments = true;

% check: if px lies inside both segments
if (((mua < 0) || (mua > 1)) || ((mub < 0) || (mub > 1)))
    check_inside_segments = false;
end

% check: if px coincides with the points
if ((norm(p1-px) < 1e-8) || (norm(p2-px) < 1e-8))
    check_inside_segments = false;
end
if ((norm(p3-px) < 1e-8) || (norm(p4-px) < 1e-8))
    check_inside_segments = false;
end


end

function [x,y,z] = get_coords_from_point(p)
x = p(1);
y = p(2);
z = p(3);
end

function [new_lattice] = split_connections(lattice)

nodes = lattice.nodes;
connectivity = lattice.connectivity;

% check if all nodes are unique
check_unique_nodes(nodes)

% empty arrays to store results
conns_deletion_index = [];
conns_to_add = [];


for c = 1:size(connectivity,1)

    % parent node ids
    n1 = connectivity(c,1);
    n2 = connectivity(c,2);

    % parent node positions
    p1 = nodes(n1,:);
    p2 = nodes(n2,:);

    % boolean array indicating collinear nodes
    idx = find_collinear_points(p1, p2, nodes);

    % set parent nodes to false
    idx(n1) = false;
    idx(n2) = false;

    % only continue if there are collinear nodes remaining
    if any(idx)

        % ids of new nodes identified
        new_node_ids = find(idx);

        % parameteric coordinates of new nodes
        alphas = zeros(size(new_node_ids));

        for i = 1:size(new_node_ids,1)
            
            % new node id
            n = new_node_ids(i);

            % new node position
            p = nodes(n, :);

            % new node parametric coordinate
            alphas(i) = compute_parameteric_point_on_line(p1, p2, p);

        end

        % sort by alphas
        [~, sort_indices] = sort(alphas);
        new_node_ids = new_node_ids(sort_indices);
        alphas = alphas(sort_indices);

        % only accept those new_nodes for which alpha in btw. 0 & 1
        accept = (alphas>0) & (alphas<1);

        % only continue if there are new nodes remaining
        if any(accept)
            
            % select acceptable node ids
            new_node_ids = new_node_ids(accept);
            
            % create new connectivity
            new_conns = [[n1; new_node_ids], [new_node_ids; n2]];

            % store new connectivity
            conns_to_add = [conns_to_add; new_conns];

            % record that this parent element now needs to be deleted
            conns_deletion_index = [conns_deletion_index, c];

        end % if any(accept)
        
    end % if any(idx)

end % loop over connectivity

% update (delete/add) connectivity
connectivity(conns_deletion_index,:) = [];
connectivity = [connectivity; conns_to_add];

% return new lattice
new_lattice = get_lattice(nodes, connectivity);

end


function [] = check_unique_nodes(nodes)
% identify and compare duplicate nodes to node with lowest id
cutoff_radius = 0.01;
idx = rgsearch(nodes, nodes, cutoff_radius);
for i=1:size(nodes,1)
    if ~all(idx{i} == i)
        error('All nodes are not unique!')
    end
end
end

function [idx] = find_collinear_points(p1, p2, p)
% Note: p is an array of size Nx3

% Ref: https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

val = cross(p-p1, p-p2);
d = sqrt(val(:,1).^2 + val(:,2).^2 + val(:,3).^2);
idx = (abs(d) < 1e-8);

end


function [coord] = find_nonzero_coord(vec)
[m,coord] = max(abs(vec));
if (abs(m)<1e-8)
    error('Zero vector provided in find_nonzero_coord()')
end
end


function [alpha] = compute_parameteric_point_on_line(p1, p2, p)
% Note: p is a vector of size 1x3
coord = find_nonzero_coord(p2-p1);
alpha = (p(coord) - p1(coord)) / (p2(coord) - p1(coord));
end

