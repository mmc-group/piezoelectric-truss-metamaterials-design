function [num]=endpointsofAonB(A,B)


[checkPtA1, onEndA1] = checkPointOnSegment([B(1,1:2);B(2,1:2)], A(1,1:2), 0);
[checkPtA2, onEndA2] = checkPointOnSegment([B(1,1:2);B(2,1:2)], A(2,1:2), 0);
[checkPtB1, onEndB1] = checkPointOnSegment([A(1,1:2);A(2,1:2)], B(1,1:2),0);
[checkPtB2, onEndB2] = checkPointOnSegment([A(1,1:2);A(2,1:2)], B(2,1:2), 0);

R=[checkPtA1;checkPtA2;checkPtB1;checkPtB2];
Q=[onEndA1;onEndA2;onEndB1;onEndB2];
Q(Q==0)=[];
R(R==0)=[];
num=length(R)-(length(Q)/2);
end