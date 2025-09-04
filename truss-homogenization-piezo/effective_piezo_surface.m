function []=effective_piezo_surface(e_mat)
% clc
% clear all
% close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%load stiffness from file:
% [Cvec]= load('Cvec.txt');
%Note: C is the 6x6 Voigt matrix
%Note Cvec is a 1x36 stiffness vector with the entries in the same order as the npj
%paper. Note that order is very important for correct plots. 
%Note: S here denotes 'Compliance', not 'Stiffness' as in the npj paper.
%Here, S is a 3x3x3x3 tensor
% e_mat=[-7.35846320530315e-10	-9.40453044143907e-13	-3.43041756466786e-10	-1.55023172911007e-11	1.31986684993260e-09	2.59194266939629e-11
% 1.62883947251996e-10	-2.58320709454446e-12	1.00226430470759e-10	3.65586382306245e-10	-2.57734423063877e-10	-1.07922012444923e-10
% -8.73171430736424e-10	-9.22074116829998e-13	-1.65343086162069e-10	-1.85900670081460e-11	1.28922148903036e-09	2.97685675204824e-11];
% e_mat=[0	0	0	0	17.0345	0
%     0	0	0	17.0345	0	0
%     -6.62281	-6.62281	23.2403	0	0	0];
% e_mat=[0	0	0	0	0	0
%     0	0	0	0	0	0
%     0	0	23.24	0	0	0];
% e_mat=[0 0 0 3.18 0 0
%     0 0 0 0 3.18 0
%     0 0 0 0 0 3.18];
e = get_tenosr(e_mat);
[Ex_pos,Ey_pos,Ez_pos,Ex_neg,Ey_neg,Ez_neg,E_pos,E_neg]=get_e_3d(e);
figure
elastSurf3(Ex_pos,Ey_pos,Ez_pos,Ex_neg,Ey_neg,Ez_neg,E_pos,E_neg)
end        
function [] = elastSurf3(Ex_pos,Ey_pos,Ez_pos,Ex_neg,Ey_neg,Ez_neg,E_pos,E_neg)
% subplot(2,1,1)
cmap_max=max(max(abs([E_pos;E_neg])));
surf(Ex_pos,Ey_pos,Ez_pos,E_pos,'EdgeColor','none');
hold on
surf(Ex_neg,Ey_neg,Ez_neg,E_neg,'EdgeColor','none');

camlight
material dull
hold on
surf( Ex_pos,Ey_pos,Ez_pos,E_pos,'FaceColor',0.5*[1,1,1],'FaceAlpha',0.05,'EdgeColor', 'none');

hold on
surf( Ex_neg,Ey_neg,Ez_neg,E_neg,'FaceColor',0.5*[1,1,1],'FaceAlpha',0.05,'EdgeColor', 'none');

axis equal
set(gca,'DataAspectRatio',[1 1 1])
bottom = [0 0 0.5];
middle=[1 1 1];
top = [0.5 0 0];
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = linspace(bottom(:,i),middle(:,i),128);
        newmap2(:,i) = linspace(middle(:,i),top(:,i),128);
    end
colormap(bluewhitered(256))
% colormap([newmap1;newmap2])
% caxis([-cmap_max,cmap_max]);

set(gcf,'Position',[360 355 553 543])
set(gca,'View',[ -29.5000   28.4000])
xlabel('$\hat e_1$')
ylabel('$\hat e_2$')
zlabel('$\hat e_3$')
set(gca,'FontSize',18)
% xlim([minE,maxE])
% ylim([minE,maxE])
% zlim([minE,maxE])
end
function [Ex_pos,Ey_pos,Ez_pos,Ex_neg,Ey_neg,Ez_neg,E_pos,E_neg] = get_e_3d(e)
n=200;
phi = linspace(0,pi,n);
theta = linspace(0,2*pi,n);
E=zeros(length(phi),length(theta));
E_pos=E;
E_neg=E;
Ex_pos=E;
Ey_pos=E;
Ez_pos=E;
Ex_neg=E;
Ey_neg=E;
Ez_neg=E;
for p = 1:length(phi)
   for t = 1:length(theta)
       [Ex_pos(p,t),Ey_pos(p,t),Ez_pos(p,t),Ex_neg(p,t),Ey_neg(p,t),Ez_neg(p,t),E_pos(p,t),E_neg(p,t)] = getE(e,phi(p),theta(t));
   end
end
% Ez=Ez./23.240;
end
function [Ex_pos,Ey_pos,Ez_pos,Ex_neg,Ey_neg,Ez_neg,E_pos,E_neg] = getE(e,phi,theta)
d=[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)];

E=0;
for i=1:3
       for j=1:3
           for k=1:3
                 E = E+e(i,j,k)*d(i)*d(j)*d(k);
           end
       end
end
if E>=0
    E_pos=E;
    E_neg=0;
Ex_pos = E*d(1);
Ey_pos = E*d(2);
Ez_pos = E*d(3);
Ex_neg = 0;
Ey_neg = 0;
Ez_neg = 0;
else
    E_neg=E;
    E=-E;
    E_pos=0;
Ex_pos = 0;
Ey_pos = 0;
Ez_pos = 0;
Ex_neg = E*d(1);
Ey_neg = E*d(2);
Ez_neg = E*d(3);
end
end

function [ee] = get_tenosr(e)
ee=zeros([3,3,3]);
for i=1:3
    for j=1:3
        for k=1:3

    if j==1 && k==1
        jk=1;
        factor=1;
    elseif j==2 && k==2
        jk=2;
        factor=1;
    elseif j==3 && k==3
        jk=3;
        factor=1;
    elseif (j==2 && k==3) || (k==2 && j==3)
        jk=4;
%         factor=1/2;
    elseif (j==1 && k==3) || (k==1 && j==3)
        jk=5;
%         factor=1/2;
    elseif (j==1 && k==2) || (k==1 && j==2)
        jk=6;
%         factor=1/2;
    end

   ee(i,j,k)=factor*e(i,jk);
         end
     end
 end        

end
