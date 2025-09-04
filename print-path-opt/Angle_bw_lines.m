function angle_AB=Angle_bw_lines(A,B)
% % Angle between two lines A and B in XY plane.
% % Endpoints of A=[x1,y1;x2,y2] ; B=[x3,y3;x4,y4]

Vec_A=(A(1,:)-A(2,:));
Vec_B=(B(1,:)-B(2,:));
costheta=dot(Vec_A,Vec_B)/(norm(Vec_A)*norm(Vec_B));
if abs(costheta-1)<1e-8
    angle_AB=0;
else
    angle_AB=acos(costheta)*180/pi;
end
end