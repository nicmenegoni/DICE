function CWxyz_range = CWxyzRange(CWxyz, CWNxyz, CWr_range, minCW_dist, maxCW_dist)
%% INFROMATION
% this function calculates the range of the space in which several circular
% window are used.
%
% Input parameters:
%   CWxyz = center of the CW plane;
%   CWNxyz = nornal vector of the circular window;
%   CWr_range = distance between the circular windows;
%   minCW_dist = minimum distance (negative) between the selected point
%   cloud and fitted circular window;
%   axCW_dist = maximum distance (positive) between the selected point
%   cloud and fitted circular window.
%
%Output parameters:
% CWxyz_range = matrix Nx3 containing all the coordinates of the centers of
%               the several circular windows.

posrange = 0 : CWr_range : maxCW_dist;
negrange = 0 : -CWr_range : minCW_dist;

if (maxCW_dist - posrange(end))> (CWr_range/2)
    posrange = [posrange, maxCW_dist];
end
if (minCW_dist - negrange(end))> (CWr_range/2)
    negrange = [negrange, minCW_dist];
end
negrange =fliplr(negrange);

range = [negrange(1:end-1),posrange]';
for rng = 1 : length(range)
CWxyz_range(rng,:) = CWxyz - range(rng) * CWNxyz;
end
end
