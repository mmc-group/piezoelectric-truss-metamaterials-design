function [Enodes,vec_map]=edge_nodes(L1,L2,L3,Coords,cnodes,ndof,vec_map)

% % Edge Numbering  (in a refernce cuboid)
% E1 : front-bottom
% E2 : back-bottom
% E3 : back-top
% E4 : front-top
% E5 : bottom-left
% E6 : bttom-right
% E7 : top-right
% E8 : top-left
% E9 : front-left
% E10: front-right
% E11: back-right
% E12: back-left


% % searching for nodes on the edges parallel to L1

% Reference Edge E1 or E5 or E9, depending on L1 in this function
Enodes1=[];
Enodes_ind=0;

for i=1:size(Coords,1)
vecA = L1/norm(L1);
vecB = Coords(i,1:3)/norm(Coords(i,1:3));
delta = abs(dot(vecA, vecB)) - 1.0;

%         angle = acos( dot(L1,Coords(i,1:3)/(norm(L1)*norm(Coords(i,1:3)))));
    if isnan(delta)
        delta=0;
    end
if abs(delta)<1e-12
    Enodes_ind=Enodes_ind+1;
    Enodes1(Enodes_ind,1)=i;
else
end
end
% Enodes1=setdiff(Enodes1,cnodes);

% % Other parallel edges

% first parallel edge
Enodes2=[];
for j=1:length(Enodes1)

    Enode_loc=Enodes1(j);

    Edge2_coords=Coords(Enode_loc,1:3)+L2;

    d2=sum((Coords(:,1:3)-Edge2_coords).^2,2);

    [val,pos]=min(d2);
    if val<1e-12
    Enodes2(j,1)=pos;
    else
    disp('error:no node found on the oppsoite edge')
    Enodes2(j,1)=0;
%     return
    end
end

% second parallel edge
Enodes3=[];
for k=1:length(Enodes1)

    Enode_loc=Enodes1(k);

    Edge3_coords=Coords(Enode_loc,1:3)+L3;

    d2=sum((Coords(:,1:3)-Edge3_coords).^2,2);

    [val,pos]=min(d2);
    if val<1e-12
    Enodes3(k,1)=pos;
    else
    disp('error:no node found on the oppsoite edge')
    Enodes3(k,1)=0;
%     return
    end
end


% third parallel edge
Enodes4=[];
for l=1:length(Enodes1)

    Enode_loc=Enodes1(l);

    Edge4_coords=Coords(Enode_loc,1:3)+L2+L3;

    d2=sum((Coords(:,1:3)-Edge4_coords).^2,2);

    [val,pos]=min(d2);
    if val<1e-12
    Enodes4(l,1)=pos;
    else
    disp('error:no node found on the oppsoite edge')
    Enodes4(l,1)=0;
%     return    
    end
end

Enodes1=setdiff(Enodes1,cnodes);
Enodes2=setdiff(Enodes2,cnodes);
Enodes3=setdiff(Enodes3,cnodes);
Enodes4=setdiff(Enodes4,cnodes);


Enodes=[Enodes1; Enodes2; Enodes3; Enodes4];

if isempty (Enodes)~=1


for iE1=1:ndof
    for E=1:size(Enodes1,1)    
        DD1=dof(Enodes1(E),iE1,ndof);   
        val1=dof(Enodes1(E),iE1,ndof);         % we use Enodes1 here because first edge is the reference edge
        vec_map=fill_vec_map(DD1,val1,vec_map);

if isempty (Enodes2)~=1
    if Enodes2(E)~=0
        DD2=dof(Enodes2(E),iE1,ndof);
        val2=dof(Enodes1(E),iE1,ndof);         % we use Enodes1 here because first edge is the reference edge
        vec_map=fill_vec_map(DD2,val2,vec_map);
    end
end

if isempty (Enodes3)~=1
    if Enodes3(E)~=0
        DD3=dof(Enodes3(E),iE1,ndof);
        val3=dof(Enodes1(E),iE1,ndof);         % we use Enodes1 here because first edge is the reference edge
        vec_map=fill_vec_map(DD3,val3,vec_map);
    end
end

if isempty (Enodes4)~=1
    if Enodes4(E)~=0
        DD4=dof(Enodes4(E),iE1,ndof);
        val4=dof(Enodes1(E),iE1,ndof);         % we use Enodes1 here because first edge is the reference edge
        vec_map=fill_vec_map(DD4,val4,vec_map);
    end
end
    end
end
else
end



end
