function [current_group_split]=split_current_group(current_group, edges_global,coords)
current_group_connectivity=edges_global(current_group,:);
current_group_nodes_indices=unique(current_group_connectivity(:));
% current_group_nodes=coords(current_group_nodes_indices,:);
% Create an adjacency matrix representing the truss connectivity
numNodes = size(coords, 1);
adjacencyMatrix = zeros(numNodes);
for i = 1:size(current_group_connectivity, 1)
    node1 = current_group_connectivity(i, 1);
    node2 = current_group_connectivity(i, 2);
    adjacencyMatrix(node1, node2) = 1;
    adjacencyMatrix(node2, node1) = 1;
end
% % reduced adjacency matrix
adjacencyMatrix1=adjacencyMatrix(current_group_nodes_indices,current_group_nodes_indices);
% Find connected components using graph analysis on edges
G = graph(adjacencyMatrix1);
bins = conncomp(G);
bins_full=zeros(1,size(coords,1));
bins_full(current_group_nodes_indices)=bins;
numGroups = max(bins);
groups = cell(1, numGroups);
group_edges = cell(1, numGroups);

for i= 1:numGroups
    edge_indices = find(bins_full(current_group_connectivity(:, 1)) == i...
        & bins_full(current_group_connectivity(:, 2)) == i);
    group_edges{i} = current_group(edge_indices);
    node_indices=find(bins_full==i);
    groups{i}.nodes = coords(node_indices, :);
    groups{i}.connectivity = current_group_connectivity(edge_indices, :);
end

% Rearrange groups based on average z-coordinate
avg_z_coordinates = zeros(numGroups, 1);
for i = 1:numGroups
    node_indices = find(bins_full==i);
    avg_z_coordinates(i) = mean(coords(node_indices, 3));
end

% Sort groups based on average z-coordinate
[~, sorted_indices] = sort(avg_z_coordinates,'ascend');

% Rearrange groups according to sorted indices
groups = groups(sorted_indices);
group_edges = group_edges(sorted_indices);

current_group_split=group_edges;


end