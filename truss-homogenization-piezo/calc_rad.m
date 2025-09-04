function [radius]=calc_rad(rho,vol,Connectivity,Coords,ndof,Nelem,Fnodes)
for e=1:Nelem
% elementDof: element degrees of freedom (Dof)
  index=Connectivity(e,:);       
  elementDof=[ndof*index(1)-6 ndof*index(1)-5 ndof*index(1)-4 ndof*index(1)-3 ... 
              ndof*index(1)-2 ndof*index(1)-1 ndof*index(1)...
              ndof*index(2)-6 ndof*index(2)-5 ndof*index(2)-4 ndof*index(2)-3 ...
              ndof*index(2)-2 ndof*index(2)-1 ndof*index(2)] ;

  x1=Coords(index(1),1);
  y1=Coords(index(1),2);
  z1=Coords(index(1),3); 
  x2=Coords(index(2),1);  
  y2=Coords(index(2),2);
  z2=Coords(index(2),3);
L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +...
         (z2-z1)*(z2-z1));
fac=1;
% % Checking if the element is an edge or face beam
f_count=0;
for Fi=1:size(Fnodes,2)
    if any(Fnodes(:,Fi)==index(1)) && any(Fnodes(:,Fi)==index(2))
    f_count=f_count+1;
    else
    end
end

if f_count==0
    fac=1;
elseif f_count==1
    fac=0.5;
elseif f_count==2
    fac=0.25;
else
    fprintf('error: element %d lies on %d faces.\n',e,f_count)
    return
end
% %   calculating the effective length of the element
L_eff=L*fac;
L_eff_arr(e)=L_eff;
end

L_eff_total=sum(L_eff_arr);
radius=sqrt(rho*vol/(pi*L_eff_total));

end