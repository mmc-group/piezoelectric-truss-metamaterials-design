clc;clear all; close all

%%% ======================================================================
% Print path planning for direct-ink-writing of truss lattices
% by Saurav Sharma, TU Delft
% Purpose: generating g-code for extrusion based printing of truss lattices
% by minimizing travel motions and ensuring no collisions
% based on: Weeks et al. (doi.org/10.1002/adma.202206958)
%%%========================================================================

% % % load coordinates of single unit cells
% unitcellname='Cube_UC';
load Unit_Cell_Coordinates/Octet_UC.mat
%%% Change the four lines at the end to set the names of files to be saved

Rot_Lat='N';
lattice.V=size(lattice.coords,1);
%%% printing nozzle diameter
d_nozzle=0.05;
%%%number of tessellations to print
tess=[1,1,1];
%%% tesselating unnit cells in X,Y,Z direcitons 'tess' times
[lattice1] = tessellate(lattice.coords, lattice.edges, tess, lattice.vecs,Rot_Lat);
lattice.coords=lattice1.nodes;
lattice.edges=lattice1.connectivity;
lattice.V=size(lattice.coords,1);
%%% for plotting tessellated lattice
%  lattice_plot (lattice.coords,lattice.edges)
%%% Determine if edge needs to be directed, depending upon its vertical inclination 
[edge_dir]=directed_edges(lattice);
tic
%%% making groups of paths that can be printed continuously
[groups]=edge_validation(lattice.coords,lattice.edges, d_nozzle);
toc
for g=1:size(groups,2)
% % out-degrees of all the vertices
edge_ids=groups{g};
graph.edges_global=lattice.edges(edge_ids,:);
vertices_ids=unique(graph.edges_global);
graph.coords=lattice.coords(vertices_ids,:);
graph=RearrangeEdges(graph,vertices_ids);
graph.V=length(vertices_ids);
for i = 1:size (graph.coords,1)
   degree(i) = sum(graph.edges==i,'all');
end
% % vertices with odd degrees
[odd_vert]=find(mod(degree,2)~=0);

% % find starting node
if isempty(odd_vert)
    start_vertex=1;
else
    for j=1:length(odd_vert)
    [is_bridge]=bridge_check(graph.edges,odd_vert(j));
        if is_bridge==0
        start_vertex=odd_vert(j);
        break
        end
    end
end

% % Adjacency table

[graph.adj]=adjacency_table(graph);
% graph.adj{4}=[1;5;3;2];
% graph.adj{5}=[4;3];
% % generating Eulerian path
graph.path=[];
for i=1:(graph.V)
    u=1;
if mod(length(graph.adj{i}),2)~=0
   u=i;
%    break
end
    
    [graph]=GenEulerPath(graph,u);
end

for edge=1:(size(graph.path,1))

    a=graph.path(edge,1);
    b=graph.path(edge,2);
    global_path{g}(edge,1)=vertices_ids(a);
    global_path{g}(edge,2)=vertices_ids(b);
end

end

% % % plotting the optimized printing path
% path_plot (lattice.coords,global_path)

path=[];
for path_ind = 1:length(global_path)
aa=global_path{path_ind};
aa=aa';
aa=aa(:);
path=[path;aa];
end
%%% write lattice coordinates and path to xlsx file
% % %  useful to save to csv for large tesselations, and then write g-code
% to avoid running algorithm for different printing parameters in gcode
% file_name='outputs/Octet_UC.xlsx';
% writematrix(lattice.coords,file_name,'Sheet','Coords')
% writematrix(path,file_name,'Sheet','Path')
% % % writing .txt into gcode format
gcode_write(path,lattice.coords)