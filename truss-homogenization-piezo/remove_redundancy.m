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