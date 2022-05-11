function StereogramEig(dipdir, dip, nplane, idx3)
% StereogamEig.m takes into account clusters and uses colour coding to
% distinguis them.
%
% Note: Original script to draw Schmidt's net is by Gerry Middleton,
% November 1995.
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland


figure('name', 'Stereogram clusters and mean orientations', 'NumberTitle', 'off')
schmidt
% Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
theta=(90-dipdir); % Poles are at the opposite direction to dip direction
r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip

% Coordinates on the strereographic projection
m=size(nplane);
for n=1:m;
xp(n) = r(n)*cosd(theta(n));
yp(n) = r(n)*sind(theta(n));

% Stereographic projection with colourcodes according to index idx3
if idx3(n)==1 
    plot(xp(n),yp(n),'sm')
elseif idx3(n)==2
    plot(xp(n),yp(n),'*r')
elseif idx3(n)==3
    plot(xp(n),yp(n),'oy')
elseif idx3(n)==4
    plot(xp(n),yp(n),'+b')
elseif idx3(n)==5
    plot(xp(n),yp(n),'^c')
elseif idx3(n)==6
    plot(xp(n),yp(n),'dg')
else 
    plot(xp(n),yp(n),'vk')
end
end
title('Stereogram Clusters and mean attitudes')
end
