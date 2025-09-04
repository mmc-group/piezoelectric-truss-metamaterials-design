close all
clear all; clc;
tic
warning ('off','all')
addpath('C:\Users\pthakolkaran\Documents\PhD\30_public_repos\npy-matlab\npy-matlab')
objective = '+e31-e33';
best = -1;
e_dataset = readtable('../generative-design/data/stiffness-vec-piezo.csv');
best_seed = 0;

x_data = readNPY(strcat('C:\Users\pthakolkaran\Documents\PhD\03_Projects\smart-trusses\generative-design\invOpt\2500\train_ptb.npy'));
adj_data = readNPY(strcat('C:\Users\pthakolkaran\Documents\PhD\03_Projects\smart-trusses\generative-design\invOpt\2500\adj_data.npy'));

e_data_fem = [];


for i=500616:500616

[lattice_temp,L1,L2,L3,E,nu,mu33,e_base] = assembleUC(x_data(i,:),adj_data(i,:));
lattice_temp.nodes=ScaleToUnit(lattice_temp.nodes,lattice_temp.connectivity);
[~,~,e,~,~]=homogenization(0.1,lattice_temp.nodes,lattice_temp.connectivity,L1,L2,L3,E,nu,mu33,e_base);
e_data_fem = [e_data_fem ; e(1,5), e(2,4), e(3,1), e(3,2), e(3,3)];

if abs(e(3,3) - e_dataset{i,5}) > 1e-10
    disp('wait')
end
% truss_plotting (lattice_temp.nodes,lattice_temp.connectivity)
end

%% ---------------------- FUNCTIONS --------------------------------------------


function [lattice_temp,L1,L2,L3,E,nu,mu33,e_base] = assembleUC(offsets,adj_triu_vals)
    nodesInit = [ [0.,0.,0.]
                [0.,0.5,0.]
                [0.,1.,0.]
                [0.5,1.,0.]
                [0.5,0.5,0.]
                [0.5,0.,0.]
                [1.,0.,0.]
                [1.,0.5,0.]
                [1.,1.,0.]
                [0.,0.,0.5]
                [0.,0.5,0.5]
                [0.,1.,0.5]
                [0.5,1.,0.5]
                [0.5,0.5,0.5]
                [0.5,0.,0.5]
                [1.,0.,0.5]
                [1.,0.5,0.5]
                [1.,1.,0.5]
                [0.,0.,1.]
                [0.,0.5,1.]
                [0.,1.,1.]
                [0.5,1.,1.]
                [0.5,0.5,1.]
                [0.5,0.,1.]
                [1.,0.,1.]
                [1.,0.5,1.]
                [1.,1.,1.]  ];

    ptb_mask = load('ptb_mask.csv');

    E=60.606e9; % Young's Modulus
    nu=0.3;  %Poisson's Ratio
    
    mu33=1433.6*8.854e-12;    %permittivity
    % piezoelectric coefficients' matrix (stress-charge form)
    e_base=[0	0	0	0	17.0345	0
        0	0	0	17.0345	0	0
        -6.62281	-6.62281	23.2403	0	0	0];
    
    
    L1=[1,0,0];
    L2=[0,1,0];
    L3=[0,0,1];
    
%     adj = round(upperTriangularToFull(adj_triu_vals),0);
    adj = round(upperTriangularToFull(adj_triu_vals),0);
    k = 1;
    ptb = zeros(27,3);
    for i = 1:size(ptb,1)
        for j = 1:size(ptb,2)
            if ptb_mask(i,j) == 1
                ptb(i,j) = offsets(k);
                k=k+1;
            end
        end
    end
    nodes = ptb2nodes(ptb);

    % % Mirror about three planes
    
    [all_nodes,adj]=mirrortruss(adj,nodes);
    
    % % % convert adjacency to connectivity table
    [connectivity]=adjtoconnect(all_nodes,adj);
    

    % % % removing duplicate nodes and edges (needed for homogenization,
    % can be omitted if only visualizing the lattices
    clear lattice_temp
    lattice_temp.nodes=all_nodes;
    lattice_temp.connectivity=connectivity;
    lattice_temp=remove_redundancy(lattice_temp);
    lattice_temp=process_lattice(lattice_temp,1);
    lattice_temp=remove_redundancy(lattice_temp);

end


function full_matrix = upperTriangularToFull(upper_triangular_values)
    n = sqrt(2 * length(upper_triangular_values) + 0.25) - 0.5; % Calculate dimension from number of elements

    full_matrix = eye(n); % Initialize full matrix with zeros

    % Fill upper triangular part of the matrix
    idx = 1;
    for i = 1:n
        for j = i:n
            full_matrix(i, j) = upper_triangular_values(idx);
            full_matrix(j, i) = upper_triangular_values(idx); % Symmetrically fill lower triangular part
            idx = idx + 1;
        end
    end
end


function result = bary_coor1D(r)
    result = r + 0.5;
end

function [ms, ns] = bary_coor2D(s, t)
    m1 = 1; n1 = 1;
    m2 = 1; n2 = 0;
    m3 = 0; n3 = 0;
    m4 = 0; n4 = 1;

    N1 = 1/4 * (1-s) * (1-t);
    N2 = 1/4 * (1-s) * (1+t);
    N3 = 1/4 * (1+s) * (1+t);
    N4 = 1/4 * (1+s) * (1-t);

    ms = N1*m1 + N2*m2 + N3*m3 + N4*m4;
    ns = N1*n1 + N2*n2 + N3*n3 + N4*n4;
end

function [xs, ys, zs] = bary_coor3D(s, t, q)
    x1 = 1; y1 = 1; z1 = 1;
    x2 = 1; y2 = 1; z2 = 0;
    x3 = 1; y3 = 0; z3 = 1;
    x4 = 1; y4 = 0; z4 = 0;
    x5 = 0; y5 = 1; z5 = 1;
    x6 = 0; y6 = 1; z6 = 0;
    x7 = 0; y7 = 0; z7 = 1;
    x8 = 0; y8 = 0; z8 = 0;

    N1 = 1/8 * (1-s) * (1-t) * (1-q);
    N2 = 1/8 * (1-s) * (1-t) * (1+q);
    N3 = 1/8 * (1-s) * (1+t) * (1-q);
    N4 = 1/8 * (1-s) * (1+t) * (1+q);
    N5 = 1/8 * (1+s) * (1-t) * (1-q);
    N6 = 1/8 * (1+s) * (1-t) * (1+q);
    N7 = 1/8 * (1+s) * (1+t) * (1-q);
    N8 = 1/8 * (1+s) * (1+t) * (1+q);

    xs = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8;
    ys = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8;
    zs = N1*z1 + N2*z2 + N3*z3 + N4*z4 + N5*z5 + N6*z6 + N7*z7 + N8*z8;
end

function nodes_n = ptb2nodes(new_ptb)
        x_edge_nodes = [4, 6, 22, 24];
        y_edge_nodes = [2, 8, 20, 26];
        z_edge_nodes = [10, 12, 16, 18];
    
        xy_face_nodes = [5, 23];
        xz_face_nodes = [13, 15];
        yz_face_nodes = [11, 17];
    
        body_nodes = [14];
        nodesInit = [ [0.,0.,0.]
        [0.,0.5,0.]
        [0.,1.,0.]
        [0.5,1.,0.]
        [0.5,0.5,0.]
        [0.5,0.,0.]
        [1.,0.,0.]
        [1.,0.5,0.]
        [1.,1.,0.]
        [0.,0.,0.5]
        [0.,0.5,0.5]
        [0.,1.,0.5]
        [0.5,1.,0.5]
        [0.5,0.5,0.5]
        [0.5,0.,0.5]
        [1.,0.,0.5]
        [1.,0.5,0.5]
        [1.,1.,0.5]
        [0.,0.,1.]
        [0.,0.5,1.]
        [0.,1.,1.]
        [0.5,1.,1.]
        [0.5,0.5,1.]
        [0.5,0.,1.]
        [1.,0.,1.]
        [1.,0.5,1.]
        [1.,1.,1.]];    
    
    nodes_n = nodesInit;

    for n = x_edge_nodes
        nodes_n(n,1) = bary_coor1D(new_ptb(n, 1));
    end
    for n = y_edge_nodes
        nodes_n(n,2) = bary_coor1D(new_ptb(n, 2));
    end
    for n = z_edge_nodes
        nodes_n(n,3) = bary_coor1D(new_ptb(n, 3));
    end

    for n = xy_face_nodes
        s = new_ptb(n,1);
        t = new_ptb(n,2);
        [nodes_n(n,1), nodes_n(n,2)] = bary_coor2D(s, t);
    end

    for n = yz_face_nodes
        s = new_ptb(n,2);
        t = new_ptb(n,3);
        [nodes_n(n,2), nodes_n(n,3)] = bary_coor2D(s, t);
    end

    for n = xz_face_nodes
        s = new_ptb(n,1);
        t = new_ptb(n,3);
        [nodes_n(n,1), nodes_n(n,3)] = bary_coor2D(s, t);
    end

    for n = body_nodes
        s = new_ptb(n,1);
        t = new_ptb(n,2);
        q = new_ptb(n,3);
        [nodes_n(n,1), nodes_n(n,2), nodes_n(n,3)] = bary_coor3D(s, t, q);
    end
end



function [mat]=SparsetoDesne(filename)
load (filename);

% Convert Python's 0-based indexing to MATLAB's 1-based indexing
indices = indices + 1;
indptr = indptr + 1;

% Determine the number of rows and columns from the indices array
numCols = max(indices); % In MATLAB, indices should represent the last index, not the size
numRows = length(indptr) - 1;

% Preallocate a dense matrix in MATLAB
denseMatrix = zeros(numRows, numCols);

% Reconstruct the dense matrix from the CSR format
for i = 1:numRows
    rowStart = indptr(i);
    rowEnd = indptr(i + 1) - 1;
    indicesForThisRow = indices(rowStart:rowEnd);
    dataForThisRow = data(rowStart:rowEnd);
    denseMatrix(i, indicesForThisRow) = dataForThisRow;
end
mat=denseMatrix;

end
function [all_nodes,adj]=mirrortruss(adj,nodes)

% Mirror about XY plane
mir_nodes=nodes;
mir_nodes(:,3)=-mir_nodes(:,3);
all_nodes=[nodes;mir_nodes];
adj=[adj,adj*0;adj*0,adj];

% Mirror about XZ plane
mir_nodes=all_nodes;
mir_nodes(:,2)=-mir_nodes(:,2);
all_nodes=[all_nodes;mir_nodes];
adj=[adj,adj*0;adj*0,adj];

% Mirror about YZ plane

mir_nodes=all_nodes;
mir_nodes(:,1)=-mir_nodes(:,1);
all_nodes=[all_nodes;mir_nodes];
adj=[adj,adj*0;adj*0,adj];
end


function [connectivity]=adjtoconnect(nodeCoordinates,adjacencyTable)

% Plot edges (truss members)
id_con=0;
for i = 1:size(adjacencyTable, 1)
    for j = i+1:size(adjacencyTable, 2)
        if adjacencyTable(i, j) == 1
            id_con=id_con+1;
            connectivity(id_con,:)=[i,j];
            % Plot a line between connected nodes
            % plot3([nodeCoordinates(i,1), nodeCoordinates(j,1)], ...
            %       [nodeCoordinates(i,2), nodeCoordinates(j,2)], ...
            %       [nodeCoordinates(i,3), nodeCoordinates(j,3)], 'b');
            % hold on
        end
    end
end
end

function [lattice] = get_lattice(nodes, connectivity)
lattice.nodes = nodes;
lattice.connectivity = connectivity;
end


function [new_lattice] = remove_redundancy(lattice)

nodes = lattice.nodes;
connectivity = lattice.connectivity;

cutoff_radius = 0.001;

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

function [lattice] = process_lattice(lattice,latt_i)

lattice = add_intersectional_nodes(lattice,latt_i);
% Do NOT run remove_redundancy() here
lattice = split_connections(lattice);

end

function [new_lattice] = add_intersectional_nodes(lattice,latt_i)
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
            [px, check_inside_segments] = intersection_point(p1, p2, p3, p4,latt_i);

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

function [px, check_inside_segments] = intersection_point(p1, p2, p3, p4,latt_i)

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
if abs(det(A))<1e-15
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
if (delta > 1e-5)
    % disp('Incorrect solution found for intersection point!')
    disp(['lattice to check' num2str(latt_i)])
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
cutoff_radius = 0.001;
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