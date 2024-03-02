function [d] = plane2point(PLANE_N, PLANE_xyz, POINT_xyz)
%% INFORMATION
% Input parameters:
%   PLANE_N = vector 1x3 representing the normal vector of the plane;
%   PLANE_XYZ = vector 1x3 containing the x, y and z coordinates of the 
%               center of the plane;
%   POINT_XYZ = vector 1x3 containing the x, y and z coordinates of the
%               point outside the plane.
% Output parameter:
%   d = distance from the plane and the point ouside the plane along the
%       normal vector of the plane.
%% SCRIPT:
% Define a unit normal vector (n) fot the normal of the plane
n = PLANE_N/norm(PLANE_N);
% define the vector(v) that join the center of the plane with the selected
% point outside the plane
%%%v = bsxfun(@minus, PLANE_xyz, POINT_xyz);
v = (POINT_xyz-PLANE_xyz);
%%%v = [PLANE_xyz(1) - POINT_xyz(1),PLANE_xyz(2) - POINT_xyz(2),PLANE_xyz(3) - POINT_xyz(3)]
% calculating the distance as the dot product of 'n' and 'v'
d = dot(n,v);
end