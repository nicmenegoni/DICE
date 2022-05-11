function Stereogram(dipdir, dip, nplane, idx3, xp, yp)
% Stereogram.m takes into account clusters and uses colour coding to
% distinguis them
%
% Note: Original script to draw Schmidt's net is by Gerry Middleton,
% November 1995.
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland

close(findobj('type','figure','name','Stereogram coloured'))
figure('name', 'Stereogram coloured', 'NumberTitle', 'off')
schmidt
% Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
theta=(90-dipdir); % Poles are at the opposite direction to dip direction
r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip

% Coordinates on the strereographic projection
m=nplane;;
for i=1:m;
xp(i) = r(i)*cosd(theta(i));
yp(i) = r(i)*sind(theta(i));

% Stereographic projection with colourcodes according to index idx3
if idx3(i)==1 
    plot(xp(i),yp(i),'sm')
elseif idx3(i)==2
    plot(xp(i),yp(i),'*r')
elseif idx3(i)==3
    plot(xp(i),yp(i),'oy')
elseif idx3(i)==4
    plot(xp(i),yp(i),'+b')
elseif idx3(i)==5
    plot(xp(i),yp(i),'^c')
elseif idx3(i)==6
    plot(xp(i),yp(i),'dg')
else 
    plot(xp(i),yp(i),'vk')
end
end
title('Stereogram Clusters')
end
