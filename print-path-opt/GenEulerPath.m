function [graph]=GenEulerPath(graph,u)
for i=1:length (graph.adj{u})
    ad=0;
    for k=1:length(graph.adj)
        ad=ad+sum(graph.adj{k});
    end
    if ad==0
        return
    end
    if i>length (graph.adj{u})
        return
    else
        v=graph.adj{u}(i);
    end
    %     intersect=edge_validation(graph,u,v)
    [valid]=is_valid(u,v,graph);
    if valid==true

        graph.path(end+1,:)=[u v];
        [graph]=rmvedge(graph,u,v);
        [graph]=GenEulerPath(graph,v);
    end

    if i==length (graph.adj{u}) && i~=1
        graph.path(end+1,:)=[u v];
        [graph]=rmvedge(graph,u,v);

    end


end

end


function [valid]=is_valid(u,v,graph)

adj_u=cell2mat(graph.adj(u));
if length(adj_u)==1
    valid = true;
else
    visited=false*ones(1,graph.V);
    %         adj_u=(graph.adj{u});
    [count1,~]=DFScount(u,visited,graph);
    [graph]=rmvedge(graph,u,v);
    visited=false*ones(1,graph.V);
    [count2,~]=DFScount(u,visited,graph);
    graph=add_edge(u,v,graph);
    if count1>count2
        valid=false;
    else
        valid=true;

    end
end
end

function [count,visited] = DFScount(v,visited,graph)
count=1;
visited(v)=true;
for i=[graph.adj{v}]'
    if visited(i)==false
        [Rcount,visited]=DFScount(i,visited,graph);
        count = count+ Rcount;

    end
end
end


function [graph]=rmvedge(graph,u,v)
adj_u=cell2mat(graph.adj(u));
adj_v=cell2mat(graph.adj(v));
for i=1:length(adj_u)
    if adj_u(i)==v
        graph.adj{u}(i)=[];
    end
end

for j=1:length(adj_v)
    if adj_v(j)==u
        graph.adj{v}(j)=[];
    end
end
end

function graph=add_edge(u,v,graph)
graph.adj{u}(end+1)=v;
graph.adj{v}(end+1)=u;
end