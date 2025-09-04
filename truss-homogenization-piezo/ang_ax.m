function [Rotax,Rotan]=ang_ax(THETA_Rot)
Rotan=THETA_Rot(1);

I=THETA_Rot(2);
J=THETA_Rot(3);
K=sqrt(1-(I^2+J^2));
K=real(K);
Rotax=[I J K];
end