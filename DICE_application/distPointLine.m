function d = distPointLine(Ppoint,Pline,Vline)
%% INFORMTATION
% This funcxtion calcualtes the distance between a line and a point in the
% 3D space. Thi script it has been developed using the solution proposed by
% Weisstein.
%References:
% Weisstein, Eric W. "Point-Line Distance--3-Dimensional." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html 
%  Input parameters:
%     Ppoint = x, y and z coordinates of the Point;
%     Pline = x, y and z coordinates of a point that below onto the line;
%     Vline = direction vector of the line;
%  Output parametrs:
%  d = distance between line and point
%% SCRIPT
p0=Ppoint;
p1=Pline;
p2=Pline+Vline;
d = abs(cross((p0-p1),(p0-p2))) / abs(p1 - p2);
d = norm(cross(Vline,(Pline-Ppoint)) /norm(Vline));
end