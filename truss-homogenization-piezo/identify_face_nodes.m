function [Fnodes_full]=identify_face_nodes(L1,L2,L3,i,j,tessaltions,Coords_tess)
if j==1
    L1=[1 0 0]*(i);
elseif j==2
    L1=[1 0 0]*(tessaltions+1);
    L2=[0 1 0]*(i);
elseif j==3
    L1=[1 0 0]*(tessaltions+1);
    L2=[0 1 0]*(tessaltions+1);
    L3=[0 0 1]*(i);
end
[cnodes]=corner_nodes_tess(L1,L2,L3,Coords_tess);
[Fnodes_full1]=face_nodes_tess(L1,L2,L3,Coords_tess,cnodes);
[Fnodes_full2]=face_nodes_tess(L2,L1,L3,Coords_tess,cnodes);
[Fnodes_full3]=face_nodes_tess(L3,L1,L2,Coords_tess,cnodes);

sizes_fnodes_full=[size(Fnodes_full1,1) size(Fnodes_full2,1) size(Fnodes_full3,1)];

Fnodes_full=zeros((max(sizes_fnodes_full)),6);
if isempty(Fnodes_full1)~=1
Fnodes_full(1:size(Fnodes_full1,1),1:2)=Fnodes_full1;
end
if isempty(Fnodes_full2)~=1
Fnodes_full(1:size(Fnodes_full2,1),3:4)=Fnodes_full2;
end
if isempty(Fnodes_full3)~=1
Fnodes_full(1:size(Fnodes_full3,1),5:6)=Fnodes_full3;
end

end
