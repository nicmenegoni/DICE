function [xp, yp]=Stereogram_colourless_Schmidt(dipdir, dip, nplane)
% Stereogram_colourless.m draws oriented data into a Schmidt's net. All the
% point are black stars.
%
% Note: Original script to draw Schmidt's net is by Gerry Middleton, November 1995. 
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland

close(findobj('type','figure','name','Stereogram'))
figure('name', 'Stereogram', 'NumberTitle', 'off')
schmidt
hold on
% Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
theta=(90-dipdir); % Poles are at the opposite direction to dip direction -> -180
r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip
m=nplane;
% Coordinates on the strereographic projection
for n=1:m;
xp(n) = r(n)*cosd(theta(n));
yp(n) = r(n)*sind(theta(n));
plot(xp(n),yp(n),'*k')
end
end