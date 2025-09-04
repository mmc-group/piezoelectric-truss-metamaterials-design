function  [KK]=...
    StiffnessMatrix_new(Sdof,ndof,Nelem,Connectivity,Coords,E1,A,Iz,Iy,Ix,...
    nu,Fnodes,mu33,e_base,e33_vec,dia1)
KK=zeros(Sdof);
Kuphi=zeros(size(Coords,1));
% computation of the system stiffness matrix
for e=1:Nelem
    dia=dia1;
    E=E1;
A=(pi/4)*dia^2; % Area of the cross section
Iy=(pi*dia^4)/64; % Moment of inertia about Y-axis
Iz=(pi*dia^4)/64; % Moment of inertia about Z-axis
Ix=(pi*dia^4)/32; % Moment of inertia about X-axis

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
% %   calculating the volume of the element
% vol=A*L*fac;
% volume_frac(e)=vol;
% EE(e,1)=fac;

% %   calculating the Shear modulus of the element
G=E/(2*(1+nu)); %Modulus of rigidity / shear modulus

k1 = E*A/L;
k2 = 12*E*Iz/(L*L*L);
k3 = 6*E*Iz/(L*L);
k4 = 4*E*Iz/L;
k5 = 2*E*Iz/L;
k6 = 12*E*Iy/(L*L*L);
k7 = 6*E*Iy/(L*L);
k8 = 4*E*Iy/L;
k9 = 2*E*Iy/L;
k10 = G*Ix/L;

a=[k1 0 0; 0 k2 0; 0 0 k6];
b=[ 0 0 0;0 0 k3; 0 -k7 0];
c=[k10 0 0;0 k8 0; 0 0 k4];
d=[-k10 0 0;0 k9 0;0 0 k5];
 
if x1 == x2 && y1 == y2
   if z2 > z1
      Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
      CXx=0;
      CYx = 0;
      CZx = 1;
      CXy = 0;
      CYy = 1;
      CZy = 0;
      CXz = -1;
      CYz = 0;
      CZz = 0;
   else
      Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
      CXx=0;
      CYx = 0;
      CZx = -1;
      CXy = 0;
      CYy = -1;
      CZy = 0;
      CXz = 1;
      CYz = 0;
      CZz = 0;

   end
else
   CXx = (x2-x1)/L;
	CYx = (y2-y1)/L;
	CZx = (z2-z1)/L;
	D = sqrt(CXx*CXx + CYx*CYx);
	CXy = -CYx/D;
	CYy = CXx/D;
	CZy = 0;
	CXz = -CXx*CZx/D;
	CYz = -CYx*CZx/D;
	CZz = D;
Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
end
[e_base_rot]=e_rotation_tensor(e_base',Lambda);
e_base_rot=e_base_rot';
e33=e_base_rot(1,1);
e31=e_base_rot(3,1);
e32=e_base_rot(3,2);
e36=e_base_rot(1,6);
% e33=e33_1*e33_vec(e);
kuu = [a b -a b;b' c b d; (-a)' b' a -b;b' d' (-b)' c];
k_phiphi=(A*mu33/L)*[1 -1
                     -1 1]; 
k_phiu=[(A*e33/L) 0 0 (Ix*e36/L) 0 0 -(A*e33/L) 0 0 -(Ix*e36/L) 0 0
        -(A*e33/L) 0 0 -(Ix*e36/L) 0 0 (A*e33/L) 0 0 (Ix*e36/L) 0 0];
k_uphi=k_phiu';
k=[kuu(1:6,1:6) k_uphi(1:6,1) kuu(1:6,7:12) k_uphi(1:6,2)
   k_phiu(1,1:6) k_phiphi(1,1) k_phiu(1,7:12) k_phiphi(1,2)
   kuu(7:12,1:6) k_uphi(7:12,1) kuu(7:12,7:12) k_uphi(7:12,2)
   k_phiu(2,1:6) k_phiphi(2,1) k_phiu(2,7:12) k_phiphi(2,2)];


R = [Lambda zeros(3,11); zeros(3) Lambda zeros(3,8); zeros(1,6) 1 zeros(1,7)
   zeros(3,7) Lambda zeros(3,4);zeros(3,10) Lambda zeros(3,1); zeros(1,13) 1];

  KK(elementDof,elementDof)=...
      KK(elementDof,elementDof)+R'*k*R*fac;
end  

end