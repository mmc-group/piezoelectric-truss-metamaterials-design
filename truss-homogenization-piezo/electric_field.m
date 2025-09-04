function [e33_vec]=electric_field(Nelem,Coords,Connectivity,VOLTAGE)

for e=1:Nelem
index=Connectivity(e,:);       
x1=Coords(index(1),1);
y1=Coords(index(1),2);
z1=Coords(index(1),3); 
x2=Coords(index(2),1);  
y2=Coords(index(2),2);
z2=Coords(index(2),3);
L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +...
         (z2-z1)*(z2-z1));
CXx = (x2-x1)/L;
CYx = (y2-y1)/L;
CZx = (z2-z1)/L;

Bphi=-[-1/L,1/L];
phi=[VOLTAGE(index(1));VOLTAGE(index(2))];
E_local(e,:)=Bphi*phi;

if E_local(e)>2e6
    e33_vec(e)=1;
elseif E_local(e)<-2e6
    e33_vec(e)=-1;
else
    e33_vec(e)=0;
end

Electric_field(e,:)=E_local(e,:)*[CXx,CYx,CZx];
avg_coords(e,1)=(x1+x2)/2;
avg_coords(e,2)=(y1+y2)/2;
avg_coords(e,3)=(z1+z2)/2;

end

Ex=Electric_field(:,1);
Ey=Electric_field(:,2);
Ez=Electric_field(:,3);
EE_net=sqrt(Ex.^2+Ey.^2+Ez.^2);
x_avg=avg_coords(:,1);
y_avg=avg_coords(:,2);
z_avg=avg_coords(:,3);
Exn=Ex./sqrt(Ex.^2+Ey.^2+Ez.^2);
Eyn=Ey./sqrt(Ex.^2+Ey.^2+Ez.^2);
Ezn=Ez./sqrt(Ex.^2+Ey.^2+Ez.^2);
scale=0.5;
% quiver3(x_avg,y_avg,z_avg,Exn,Eyn,Ezn,scale,'color','r','LineWidth',2);

end