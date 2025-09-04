function [checkPt, onEnd] = checkPointOnSegment(segment, C, plotResults)

% Check if a point C lies on a line segment
% Original author: Vipul Lugade
% Matlabgeeks.com
% 9/9/2018
% adapted by: Saurav Sharma, TU Delft
eps=1e-8;
%% initialize outputs
checkPt = false;
onEnd = false;

%% Error check
if any(size(segment) ~=2) || numel(C) ~= 2
    error('Point/segment arguments not properly defined.');
end

%% Perform cross product to see if the vectors formed by both endpoints
% and the point in question are collinear
A = segment(1,:);
B = segment(2,:);

% Check if the segment degenerates to a single point
if isequal(A, B) && isequal(A, C)
    checkPt = true;
    onEnd = true;
    return
elseif isequal(A, B) && ~isequal(A, C)
    return % Point C doesn't lie on the degenerate segment
end

% form vectors for the line segment (AB) and the point to one endpoint of
% segment
AB = B - A;
AC = C - A;
XX=cross(AB,AC);
% if not collinear then return false as point cannot be on segment
if abs(XX)<eps
    % calculate the dotproduct of (AB, AC) and (AB, AB) to see point is now
    % on the segment
    dotAB = dot(AB, AB);
    dotAC = dot(AB, AC);
    % on end points of segment
    if abs(dotAC) < eps || abs(dotAC-dotAB) <eps
        onEnd = true;
        checkPt = true;
        % on segment
    elseif dotAC > eps && dotAC < dotAB
        checkPt = true;
    end
end

%% plot the results
if plotResults
    figure; hold on
    plot(segment(:,1), segment(:,2), 'ro-');
    plot(C(1), C(2), 'bs');
    legend('segment', 'point');
    if checkPt
        text(C(1)+0.1, C(2)+0.1, 'Point on Segment');
    else
        text(C(1)+0.1, C(2)+0.1, 'Point Not on Segment');
    end
end
end
%% Cross product returning z value
function z = cross(a, b)
z = a(1)*b(2) - a(2)*b(1);
end