function [T,cnodes,Enodes,Fnodes,Fnodes_full]=PBC_mat(ndof,Sdof,L1,L2,L3,Coords)

% % Generating vector map 
% 
% each entry in vec_map denotes the dof to which
% the dof associated with the position in vec_map is coupled with 
% e.g. vec_map(3)=2 means that 3rd dof is coupled to 2nd dof

vec_map=-ones(Sdof,1);

% % Indetifying and coupling corner nodes with the reference node
[cnodes,vec_map]=corner_nodes(L1,L2,L3,Coords,ndof,vec_map);

% % Coupling Edges with the edges E1, E5, E9 (Description of edge notation can be found in the function "edge_nodes")
[Enodes1,vec_map]=edge_nodes(L1,L2,L3,Coords,cnodes,ndof,vec_map);
[Enodes2,vec_map]=edge_nodes(L2,L1,L3,Coords,cnodes,ndof,vec_map);
[Enodes3,vec_map]=edge_nodes(L3,L2,L1,Coords,cnodes,ndof,vec_map);
sizes_Enodes=[length(Enodes1) length(Enodes1) length(Enodes1)];
Enodes=zeros(max(sizes_Enodes));
Enodes(1:length(Enodes1),1)=Enodes1;
Enodes(1:length(Enodes2),2)=Enodes2;
Enodes(1:length(Enodes3),3)=Enodes3;

% % Coupling Faces with the Faces F1, F2, F3 (Description of face notation can be found in the function "face_nodes")
[Fnodes1,Fnodes_full1,vec_map]=face_nodes(L1,L2,L3,Coords,cnodes,Enodes,ndof,vec_map);
[Fnodes2,Fnodes_full2,vec_map]=face_nodes(L2,L1,L3,Coords,cnodes,Enodes,ndof,vec_map);
[Fnodes3,Fnodes_full3,vec_map]=face_nodes(L3,L2,L1,Coords,cnodes,Enodes,ndof,vec_map);

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

sizes_fnodes=[length(Fnodes1) length(Fnodes1) length(Fnodes3)];

Fnodes=zeros((max(sizes_fnodes)),3);
if isempty(Fnodes1)~=1
Fnodes(1:length(Fnodes1),1)=Fnodes1;
end
if isempty(Fnodes2)~=1
Fnodes(1:length(Fnodes2),2)=Fnodes2;
end
if isempty(Fnodes3)~=1
Fnodes(1:length(Fnodes3),3)=Fnodes3;
end

% Fnodes=[Fnodes1 Fnodes2 Fnodes3];
vec_map=fnalize_vec_map(vec_map);


% % Generating T matrix from vec_map 
m=length(unique(vec_map));
T=zeros (Sdof,m);

Q=zeros(Sdof,Sdof);

for i=1:Sdof
    Q(i,vec_map(i))=1;
end

S=sum(Q,1);

T_ind=find(S);

T=Q(:,T_ind);
end


