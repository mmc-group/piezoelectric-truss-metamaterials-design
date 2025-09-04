clc;clear all; close all
vertices=[0,0,0
1,0,0
1,1,0
0,1,0
0,0,1
1,0,1
1,1,1
0,1,1
0.5,0,0.5
0.5,1,0.5
0,0.5,0.5
1,0.5,0.5
0.5,0.5,0
0.5,0.5,1];
% edges=[1,13
% 2,13
% 3,13
% 4,13
% 5,14
% 6,14
% 7,14
% 8,14
% 1,11
% 4,11
% 5,11
% 8,11
% 2,12
% 3,12
% 6,12
% 7,12
% 1,9
% 2,9
% 5,9
% 6,9
% 3,10
% 4,10
% 7,10
% 8,10
% 9,11
% 9,12
% 9,13
% 9,14
% 10,11
% 10,12
% 10,13
% 10,14
% 11,13
% 11,14
% 12,13
% 12,14];
edges=1+[0,1
       0,2
       1,2
       2,3];
% lattice_plot (vertices,edges)
edge_table=table(edges,'VariableNames',{'EndNodes'})
g=graph(edge_table);
plot(g)
% % out-degrees of all the vertices
for i = 1:size (vertices)
   degree(i) = sum(edges==i,'all')
end
[odd_vert]=find(mod(degree,2)~=0);
% % vertices with odd degrees
% % generating Eulerian path

% % find starting node
if isempty(odd_vert)
    start_vertex=1;
else
    for j=1:length(odd_vert)
    [is_bridge]=bridge_check(edges,odd_vert(j))
        if is_bridge==0
        start_vertex=odd_vert(j);
        break
        end
    end
end

% % Fluery algorithm


count=1;
ind=1; 
while ind == 1
    Euler_path(count)=start_vertex;
    count=count+1;
    v=dfsearch(g,start_vertex)

    

end