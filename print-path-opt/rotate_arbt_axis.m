function [R]=rotate_arbt_axis(axis,theta)

x=axis(1);
y=axis(2);
z=axis(3);
c = cos(theta);
s = sin(theta);
t = 1 - c;

R =  [t*x*x + c	  t*x*y - z*s	   t*x*z + y*s
      t*x*y + z*s	  t*y*y + c	       t*y*z - x*s
      t*x*z - y*s	  t*y*z + x*s	   t*z*z + c];
end