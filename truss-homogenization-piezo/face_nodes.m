function [Fnodes,Fnodes_full,vec_map]=face_nodes(L1,L2,L3,Coords,cnodes,Enodes,ndof,vec_map)
% % Face Numbering  (in a refernce cuboid)
% F1 : Left
% F2 : Front
% F3 : Bottom
% F4 : right
% F5 : Back
% F6 : top

% % searching for nodes on the face F; where L1 is normal to F

Fnodes1=[];
Fnodes_ind=0;

for i=1:size(Coords,1)
     Triple_prod=dot(Coords(i,1:3),(cross(L2,L3)));
    if abs(Triple_prod)<=1e-12
    Fnodes_ind=Fnodes_ind+1;
    Fnodes1(Fnodes_ind,1)=i;
    end
end
% Fnodes1=setdiff(Fnodes1,cnodes);
% Fnodes1=setdiff(Fnodes1,Enodes);

% % Other parallel face

Fnodes2=[];

for j=1:length(Fnodes1)

    Fnode_loc=Fnodes1(j);   

    Face2_coords=Coords(Fnode_loc,1:3)+L1;

    d2=sum((Coords(:,1:3)-Face2_coords).^2,2);

    [val,pos]=min(d2);
    if val<1e-12
    Fnodes2(j,1)=pos;
    else
    disp('error:no node found on the oppsoite face')
    Fnodes2(j,1)=0;
    end
end

Fnodes_full=[Fnodes1 Fnodes2];
Fnodes1=setdiff(Fnodes1,cnodes);
Fnodes1=setdiff(Fnodes1,Enodes);
Fnodes2=setdiff(Fnodes2,cnodes);
Fnodes2=setdiff(Fnodes2,Enodes);


Fnodes=[Fnodes1; Fnodes2];

if isempty(Fnodes)~=1
for iF1=1:ndof
    for F=1:size(Fnodes1,1)    
        DD1=dof(Fnodes1(F),iF1,ndof);
        val1=dof(Fnodes1(F),iF1,ndof);         % we use Fnodes1(F) here because Fnodes1 is the reference face
        vec_map=fill_vec_map(DD1,val1,vec_map);

    if Fnodes2(F)~=0
        DD2=dof(Fnodes2(F),iF1,ndof);
        val2=dof(Fnodes1(F),iF1,ndof);         % we use Fnodes1(F) here because Fnodes1 is the reference face
        vec_map=fill_vec_map(DD2,val2,vec_map);
    end

    end
end
else
end
end