function graph=RearrangeEdges(graph,vertices_ids)
graph.edges=[];
for i=1:size(graph.edges_global,1)
    a=graph.edges_global(i,1);
    b=graph.edges_global(i,2);
    graph.edges(i,1)=find(vertices_ids==a);
    graph.edges(i,2)=find(vertices_ids==b);
end
end