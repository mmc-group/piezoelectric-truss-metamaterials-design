function [cnodes,vec_map]=corner_nodes(L1,L2,L3,Coords,ndof,vec_map)

cnodes=[];

% % Identifying the reference node
R_node=[0 0 0];
C=zeros(8,3);
% % Expected corner node coordinates
C(1,:)=R_node;
C(2,:)=R_node+L1;
C(3,:)=R_node+L2;
C(4,:)=R_node+L3;
C(5,:)=R_node+L1+L2;
C(6,:)=R_node+L1+L3;
C(7,:)=R_node+L2+L3;
C(8,:)=R_node+L1+L2+L3;

% % searching for corner nodes in the nodes list
for i=1:size(C,1)
d2=sum((Coords(:,1:3)-C(i,:)).^2,2);
[val,pos]=min(d2);
if val<1e-8
    cnodes(i,1)=pos;
else
end

end

cnodes(find(cnodes==0))=[];
% % Coupling corner nodes with the refernce (origin) node

for ic=1:ndof
    for c=1:size(cnodes,1)    
        DD=dof(cnodes(c),ic,ndof);
        val=dof(cnodes(1),ic,ndof); % we use cnodes(1) here because first node is the reference/origin node
        vec_map=fill_vec_map(DD,val,vec_map);

    end
end
end