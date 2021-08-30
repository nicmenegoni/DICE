function [A] = intCirc2Circ(c1, r1, c2, r2)
%%
% 24th April 2020 : correction of the area calculation for the second case
%% INFORMATION:
% This function calculates the area of the intersection between two
% 3D circles lying onto the same 3D plane
%
%   INPUT PARAMETERS
%   r1 and c1 are respectively the radius and the center coordinates (x, y
%   and z) of the Circle 1 and the r2 and c2 are referred to Circle 2.
%
%   CALCULATION PARAMETRS (present only inside this function)
%   d = distance between the centers of the circles
%
%   OUTPUT PARAMETERS
%   A is the area of intersection between two circles

%   Giving 2 3D discs lying onto the same 3D plane, 3 cases could happen.
%
%   -1st case, one circle (circle min) is totally included into the other
%   (circle max) and therefore the area of intersection is equal to the 
%   area of the included circle:
%
%                      d < r_max     &       d + r_min < r_max ;
%
%   -2nd case, the circles intersect each other sharing some portion of 
%   them:
%
%                  a)    d < r1      &       d + r2 > r1 ;
%                                  or
%                  b)    d > r1      &       d - r2 < r1 ;
%
%   -3rd case, the circles do not intersect each other.
%
%                      d > r1      &       d - r2 > r1 ;
%
%NB. When using this function in the P32 calculation code, the c1 and r1
%are referred to the scan circle

%% SCRIPT
% It starts calculating the distance between the center of the circles:
% d = distance between the centers of the circles
d = sqrt((c1(1)-c2(1))^2 ...
    + (c1(2)-c2(2))^2 ...
    + (c1(3)-c2(3))^2);

if d>(r1+r2) || d==(r1+r2) %3rd case, no intersection;
    A=0;
else
    if r1>r2 && ((d+r2)<r1  || (d+r2)==r1)% 1st case: circ2 in circ 1
        A = pi * r2^2;
    elseif r1<r2 && (d<r2) && ((d+r1)<r2  || (d+r1)==r2) % 1st case: circ1 in circ 2
        A = pi * r1^2;
    elseif d==0 && r1==r2% 1st case: circ1 coincides with circ 2
        A = pi * r1^2;
    else % 2nd case
        d1 = (r1^2 - r2^2 + d^2)/(2*d);
        d2 = d - d1;
        A = r2^2 * acos(d2/r2) - d1 * sqrt(r1^2 - d1^2) ...
            + r1^2 * acos(d1/r1) - d2 * sqrt(r2^2 - d2^2);
    end
end


end