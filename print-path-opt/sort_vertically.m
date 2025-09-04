function [lattice]=sort_vertically(lattice)
c=lattice.coords;
n=lattice.edges;

[z,index]=sort(c(:,3));
c=c(index,:);
for i=1:size(n,1)
    a=find(index==n(i,1));
    b=find(index==n(i,2));
    
n(i,1)=a;
n(i,2)=b;
end
lattice.coords=c;
lattice.edges=n;
end
