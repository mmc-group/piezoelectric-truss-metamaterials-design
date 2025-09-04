function [Pz]=intersection_point(A,Px,Py)
x1=A(1,1);y1=A(1,2);z1=A(1,3);x2=A(2,1);y2=A(2,2);z2=A(2,3);

l=(x2-x1);m=(y2-y1);n=(z2-z1);
Pz=(((Px-x1)*n)/l)+z1;
end