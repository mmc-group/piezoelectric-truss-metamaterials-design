function [adj]=adjacency_table(graph)

for i = 1:size(graph.coords,1)
ad=[];
    [a,b]=find(graph.edges==i);
for j=1:length(a)
if b(j)==1
    ad(j)=graph.edges(a(j),2);
elseif b(j)==2
    ad(j)=graph.edges(a(j),1);
end
end
    adj{i}=sort(ad');
end