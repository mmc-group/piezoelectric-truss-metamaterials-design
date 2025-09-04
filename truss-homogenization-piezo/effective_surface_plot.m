function []=effective_surface_plot(Cvec)
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
% Cvec=[1.27205e+011 	8.02122e+010 	8.46702e+010 	0 	0 	0 
% 8.02122e+010 	1.27205e+011 	8.46702e+010 	0 	0 	0 
% 8.46702e+010 	8.46702e+010 	1.17436e+011 	0 	0 	0 
% 0 	0 	0 	2.29885e+010 	0 	0
% 0	0	0	0	2.29885e+010	0
% 0	0	0	0	0	2.34742e+010];

S = getCompliance(Cvec);
figure
elastSurf3(S)
end        
function [] = elastSurf3(S)
n=200;
phi = linspace(0,pi,n);
theta = linspace(0,2*pi,n);
[Ex,Ey,Ez,E] = getE_3d(S,phi,theta);
Ex=Ex*1e-9;
Ey=Ey*1e-9;
Ez=Ez*1e-9;
E=E*1e-9;

surf(Ex,Ey,Ez,E,'EdgeColor','none');
minE = min(min(E));
maxE = max(max(E));
caxis([minE,maxE]);
camlight
material dull
hold on
surf( Ex,...
      Ey,...
      Ez,...
    'FaceColor',0.5*[1,1,1],'FaceAlpha',0.05,'EdgeColor', 'none');
axis equal
set(gca,'DataAspectRatio',[1 1 1])
colormap winter
set(gcf,'Position',[360 255 453 443])
set(gca,'View',[ -29.5000   28.4000])
xlabel('$\hat e_1$')
ylabel('$\hat e_2$')
zlabel('$\hat e_3$')
set(gca,'FontSize',18)
xlim([-maxE,maxE])
ylim([-maxE,maxE])
zlim([-maxE,maxE])
end
function [Ex,Ey,Ez,E] = getE_3d(S,phi,theta)
E=zeros(length(phi),length(theta));
Ex=E;
Ey=E;
Ez=E;
for p = 1:length(phi)
   for t = 1:length(theta)
       [Ex(p,t),Ey(p,t),Ez(p,t),E(p,t)] = getE(S,phi(p),theta(t));
   end
end
end
function [x,y,z,E] = getE(S,phi,theta)
d=[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)];
E = (1/doubleContraction(S,d));
x = E*d(1);
y = E*d(2);
z = E*d(3);
end
function [x] = doubleContraction(S,d)
x=0;
for i=1:3
       for j=1:3
           for k=1:3
               for l=1:3
                   x = x+S(i,j,k,l)*d(i)*d(j)*d(k)*d(l);
               end
           end
       end
end
end
   
function [S] = getCompliance(vec)
map = [...
    1,6,5;...
    6,2,4;...
    5,4,3];
Cv = reshape(vec,[6,6]);
Cv = 0.5*(Cv+Cv');
Sv = inv(Cv);
S=zeros([3,3,3,3]);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                vr = map(i,j);
                vc = map(k,l);
                
                factor=1;
                if(vr>3)
                    factor = 2*factor;
                end
                if(vc>3)
                    factor = 2*factor;
                end
                
                S(i,j,k,l) = 1/factor * Sv(vr,vc);
                    
            end
        end
    end
end
end