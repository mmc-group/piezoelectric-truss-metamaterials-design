function [edge_dir]=directed_edges(graph)
edge_dir=zeros(length(graph.edges),1);
z=[0,0,1];
% connectivity=graph.edges;
for i =1:length(graph.edges)
edge=[graph.edges(i,1),graph.edges(i,2)];
A=graph.coords(edge(1),:);
B=graph.coords(edge(2),:);
AB=(A-B);
costheta=dot(AB,z)/(norm(AB)*norm(z));
theta=acos(costheta)*180/pi;
if theta<=30
    edge_dir(i)=1;
else
    edge_dir(i)=0;
end
end
end