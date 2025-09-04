function [Relative_density,C,e,mu,d]=homogenization(rho,Coords,Connectivity,L1,L2,L3,E,nu,mu33,e_base)
% % Mesh Parameters
xx=Coords(:,1);
yy=Coords(:,2);
zz=Coords(:,3);
Nnodes=size(Coords,1);
Nelem=size(Connectivity,1);
ndof=7;
Sdof=ndof*Nnodes; 

% % Plotting the truss lattice
% truss_plotting(Coords,Connectivity)

Volume=dot(L1,cross(L2,L3));

% % % Generating Trnasformation matrix for Periodic boundary condition

[T,~,~,~,Fnodes_full]=PBC_mat(ndof,Sdof,L1,L2,L3,Coords);
[radius]=calc_rad(rho,Volume,Connectivity,Coords,ndof,Nelem,Fnodes_full);
d=2*radius;
% d=0.05*2;
A=(pi/4)*d^2; % Area of the cross section
Iy=(pi*d^4)/64; % Moment of inertia about Y-axis
Iz=(pi*d^4)/64; % Moment of inertia about Z-axis
Ix=(pi*d^4)/32; % Moment of inertia about X-axis



% % Translation of the entire unit cell
% offset=[-0.5,-0.5,-0.5];
% Coords(:,1:3)=Coords(:,1:3)+offset;



% % % Determining sign of e33 based on poling direction in Local 
% % % coordinates of individual struts
e33_vec=[];
% [e33_vec]=poling_direction(Nelem,Connectivity,Coords,A,mu33,L3);
% % % Tessellating the unit cell
% tessaltions=1;
% [displ,C11]=FEM_load_case(Coords(:,1:3),Connectivity,tessaltions,L1,L2,...
%     L3,Fnodes_full,E,d,nu,mu33,e36,e33,d,dd1);
[K]=StiffnessMatrix(Sdof,ndof,Nelem,Connectivity,...
    Coords,E,A,Iz,Iy,Ix,nu,Fnodes_full,mu33,e_base,e33_vec,d);

% [U]=FEM(K,Coords(:,1:3),Connectivity,Sdof,ndof)
% [E]=electric_field(K,Nelem,Coords,Connectivity,Sdof,xx,yy,zz)
Volume=dot(L1,cross(L2,L3));
Relative_density=rho;

% % Homogenization
[C,e,mu]=Effective_stiffness(T,K,Sdof,xx,yy,zz,Fnodes_full);
C=C/Volume;
e=e/Volume;
mu=mu/Volume;
end