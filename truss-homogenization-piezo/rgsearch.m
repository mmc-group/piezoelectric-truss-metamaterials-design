function [idx_all,dist_all]=rgsearch(X,Y,r)
n = size(Y,1);
idx_all = cell(n,1);
dist_all = cell(n,1);
for i=1:n
    [idx, dist] = rgsearch_single(Y(i,:),r,X,2);
    [~,sorted_indices] = sort(dist);
    idx = idx(sorted_indices);
    dist = dist(sorted_indices);
    idx_all{i} = idx';
    dist_all{i} = sqrt(dist');
end

end

function [idx,dist]=rgsearch_single(c,r,X,mode)
% RANGESEARCH Range Search to find all points within a range.
% [index,distance]=rangesearch(c,r,X,mode) returns all points,x of X which
% are in the range of ||x-c||<r.
% Inputs:
%    c: 1 x d the querry vector
%    r: scalar defines the range
%    X: n x d array of all search points
% mode: range mode, 1 - box range, 2 (default) - radial range 
% Outputs:
%    index: indices of points in the range.
% distance: distances between the reference point and points in the range.
%
% See Also: pdist, kdtree, knnsearch

% Version 2.0 by Yi Cao at Cranfield University on 6th April 2008
%

%Examples
%Example 1: Radial range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X);
t=0:0.02:2*pi;
x=c(1)+r*cos(t);
y=c(2)+r*sin(t);
subplot(211)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Radial range search')
%}
%Example 2: Box range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X,1);
t=0:0.02:2*pi;
x=c(1)+[-r r r -r -r];
y=c(2)+[r r -r -r r];
subplot(212)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Box range search')
%}
% Example 3: Large data set search
%{
N=250000;
d=10;
X=randn(N,d);
c=randn(1,d);
r=1;
tic
idx=rangesearch(c,r,X);
toc
%}
% The time is similar to using kdtree (FEX ID 4586).

if nargin<4
    mode=2;
end
if mode==2
    [idx,dist]=range2(c,r,X);
else
    [idx,dist]=range1(c,r,X);
end
end

function [idx,dist]=range1(c,r,X)
[nPoints,nVariables]=size(X);
s=zeros(nPoints,1);
for d=1:nVariables
    x=abs(X(:,d)-c(d));
    id=x>s;
    s(id)=x(id);
end
fidx=s<r;
idx=find(fidx);
dist=s(fidx);
end

function [idx,dist]=range2(c,r,X)
nVariables=numel(c);
r2=r*r;
s=0;
for d=1:nVariables
    s=s+(X(:,d)-c(d)).^2;
end
fidx=s<r2;
idx=find(fidx);
dist=s(fidx);
end
