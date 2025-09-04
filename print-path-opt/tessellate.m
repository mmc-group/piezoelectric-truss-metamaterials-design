function [new_lattice] = tessellate(nodes, connectivity, tess, lattice_vecs,Rot_Lat)

% nodes, connectivity, tessellation, lattice vectors 1, 2, & 3
% Note: lattice vectors must be for the ** untessellated ** cell 

if all(tess==1)
    % do nothing
else
% divide by "tess" to keep the whole strucutre unit dimensional
    vec1 = lattice_vecs{1};%/tess(1);
    vec2 = lattice_vecs{2};%/tess(2);
    vec3 = lattice_vecs{3};%/tess(3);

    n_unit = nodes;%./tess;
    c_unit = connectivity;

    nodes=[];
    connectivity=[];

    for i=1:tess(1)
        for j=1:tess(2)
            for k=1:tess(3)
                
                %order of next two lines is critical:
                connectivity = [connectivity; c_unit + size(nodes,1)];
                nodes = [nodes; (n_unit + (i-1)*vec1 + (j-1)*vec2 + (k-1)*vec3)];
            end
        end
    end

end

lattice = get_lattice(nodes, connectivity);
[new_lattice] = remove_redundancy(lattice);
%% Rotating lattice
if Rot_Lat=='Y'
axis=[1,0,0];
theta=pi/4;
[R]=rotate_arbt_axis(axis,theta);

rot_coords=R*new_lattice.nodes';
rot_coords=rot_coords';
new_lattice.nodes=rot_coords;
else
end
% lattice_print = extract_UnitCube(new_lattice);
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


function [idx_all,dist_all]=rgsearch(X,Y,r)
n = size(Y,1);
idx_all = cell(n,1);
dist_all = cell(n,1);
for i=1:n
    [idx, dist] = rgsearch_single(Y(i,:),r,X,2);
    [~,sorted_indices] = sort(dist);
    idx = idx(sorted_indices);
    dist = dist(sorted_indices);
    idx_all{i} = idx';
    dist_all{i} = sqrt(dist');
end

end

function [idx,dist]=rgsearch_single(c,r,X,mode)
% RANGESEARCH Range Search to find all points within a range.
% [index,distance]=rangesearch(c,r,X,mode) returns all points,x of X which
% are in the range of ||x-c||<r.
% Inputs:
%    c: 1 x d the querry vector
%    r: scalar defines the range
%    X: n x d array of all search points
% mode: range mode, 1 - box range, 2 (default) - radial range 
% Outputs:
%    index: indices of points in the range.
% distance: distances between the reference point and points in the range.
%
% See Also: pdist, kdtree, knnsearch

% Version 2.0 by Yi Cao at Cranfield University on 6th April 2008
%

%Examples
%Example 1: Radial range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X);
t=0:0.02:2*pi;
x=c(1)+r*cos(t);
y=c(2)+r*sin(t);
subplot(211)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Radial range search')
%}
%Example 2: Box range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X,1);
t=0:0.02:2*pi;
x=c(1)+[-r r r -r -r];
y=c(2)+[r r -r -r r];
subplot(212)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Box range search')
%}
% Example 3: Large data set search
%{
N=250000;
d=10;
X=randn(N,d);
c=randn(1,d);
r=1;
tic
idx=rangesearch(c,r,X);
toc
%}
% The time is similar to using kdtree (FEX ID 4586).

if nargin<4
    mode=2;
end
if mode==2
    [idx,dist]=range2(c,r,X);
else
    [idx,dist]=range1(c,r,X);
end
end

function [idx,dist]=range1(c,r,X)
[nPoints,nVariables]=size(X);
s=zeros(nPoints,1);
for d=1:nVariables
    x=abs(X(:,d)-c(d));
    id=x>s;
    s(id)=x(id);
end
fidx=s<r;
idx=find(fidx);
dist=s(fidx);
end

function [idx,dist]=range2(c,r,X)
nVariables=numel(c);
r2=r*r;
s=0;
for d=1:nVariables
    s=s+(X(:,d)-c(d)).^2;
end
fidx=s<r2;
idx=find(fidx);
dist=s(fidx);
end

function lattice_print = extract_UnitCube(new_lattice)
centre_x=(max(new_lattice.nodes(:,1))+min(new_lattice.nodes(:,1)))/2;
centre_y=(max(new_lattice.nodes(:,2))+min(new_lattice.nodes(:,2)))/2;
centre_z=(max(new_lattice.nodes(:,3))+min(new_lattice.nodes(:,3)))/2;

edge_length=0.5;
corner=[centre_x-edge_length/2 centre_y-edge_length/2 centre_z-edge_length/2];
xmax=centre_x+edge_length/2;
ymax=centre_y+edge_length/2;
zmax=centre_z+edge_length/2;


end
