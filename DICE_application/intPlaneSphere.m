function circle = intPlaneSphere(plane, sphere)
%INTERSECTPLANESPHERE Return intersection circle between a plane and a sphere
%%  INFORMATIONS

%   This is a modified version of the function 'intersectPlaneSphere' of the
%   functions package 'geom3d-2016' of David Legland.
%   The theta, phi and psi variables have been erased because not usefull in
%   the P32 calculation code.
%   This function returns the circle which is the intersection of the given 
%   plane and sphere. 
%
%   INPUT PARAMETERS
%   CIRC = intersectPlaneSphere(PLANE, SPHERE)
%   
%   PLANE  : [x0 y0 z0  Nx Ny Nz] (*modified by NM)
%   SPHERE : [XS YS ZS  RS]
%   CIRC   : [XC YC ZC  RC ] (*modified by NM)
%   [x0 y0 z0] is the origin of the plane, [dx1 dy1 dz1] and [dx2 dy2 dz2]
%   are two direction vectors,
%   [XS YS ZS] are coordinates of the sphere center, RS is the sphere
%   radius, 
%   [XC YC ZC] are coordinates of the circle center, RC is the radius of
%   the circle.
%   
%   See Also:
%   planes3d, spheres, circles3d, intersectLinePlane, intersectLineSphere
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%
%   modified by Niccolò Menegoni (NM)
%   Dept. Earth and Environmental Sciences, University of Pavia (Italy)

%%  HISTORY
%   27/06/2007: change output format of circle, add support for multiple
%       data
%   2011-06-21: use degrees for angles
%   2019-01-29: modify by NM the circle output, the angles are been erased
%   2020-05-20: modify by NM, simplfy calculation (elimination bsxfun)


%1) center of the spheres
center  = sphere(1,1:3); %modified by NM
Rs  = sphere(1,4); % modified by NM


%2) Projection of sphere center on plane -> gives circle center
origins(:) = plane(1,1:3); %modified by NM
normals(:) = plane(1,4:6); %Modified by NM

%2.1) difference between origins of plane and point
dp = origins - center;%Modified by NM

% 2.2) relative position of point on normal's line
t = sum(normals .* dp,2) ./ sum(normals.^2,2);%Modified by NM

%2.3) add relative difference to project point back to plane
circle0 = center + t .* normals;%Modified by NM

%4) radius of circles
d = center - circle0; % row modified by Niccolò Menegoni
d = sqrt(d(1)^2 + d(2)^2 + d(3)^2);%Modified by NM
Rc  = sqrt(Rs.*Rs - d.*d);%Modified by NM

%5) circle
circle = [circle0 Rc];%Modified by NM
