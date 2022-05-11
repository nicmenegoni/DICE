function Unitvectors_coloured(dipdir, dip, nplane, idx3)
% Unitvectors_coloured.m draws colour coded unit vectors on a stereographic sphere, in which dip 0
% degrees is at the bottom of the sphere and 90 degrees at the top.
%
% 12.10.20101 M. Markovaara-Koivisto, Aalto University School of Science
% and Technology, Finland

m=nplane;
theta=(dipdir+180); %Polevectors are at the opposite direction to the dip direction -> -180
X3=[(dip-45)/45]'; %z-coordinate on the sterepgraphic sphere
r=cos(asin(X3)); % Poles are perpendicular to the dip
for n=1:m;
X1(n) = r(n)*sind(theta(n)); % x-coordinate
X2(n) = r(n)*cosd(theta(n)); % y-coordinate
end

% Colour coding according to clustering
figure('name', 'Unitvectors coloured', 'NumberTitle', 'off')
hold on
m=nplane;
for n=1:m;
if idx3(n)==1 
    plot3(X1(n),X2(n),X3(n),'sm')
elseif idx3(n)==2
    plot3(X1(n),X2(n),X3(n),'*r')
elseif idx3(n)==3
    plot3(X1(n),X2(n),X3(n),'oy')
elseif idx3(n)==4
    plot3(X1(n),X2(n),X3(n),'+b')
elseif idx3(n)==5
    plot3(X1(n),X2(n),X3(n),'^c')
elseif idx3(n)==6
    plot3(X1(n),X2(n),X3(n),'dg')
else
    plot3(X1(n),X2(n),X3(n),'vk')
end
end
[xi,yi,zi]=sphere(36); % Sphere as the axes
mesh(xi,yi,zi, 'FaceColor', 'none', 'EdgeColor', 'k')
axis equal
grid on
args={'BackgroundColor','none','FontSize',8};
h(1)=text(-1.2, 0,'WEST', 'VerticalAlignment','middle', 'HorizontalAlignment','left', args{:});
h(2)=text( 1.1, 0,'EAST', 'VerticalAlignment','middle', 'HorizontalAlignment','left',args{:});
h(3)=text( 0,-1.1,'SOUTH','VerticalAlignment','middle', 'HorizontalAlignment','left', args{:});
h(4)=text( 0, 1.1,'NORTH','VerticalAlignment','middle', 'HorizontalAlignment','left', args{:});
v(5)=text( 0, 0, 1.1,'Dip 90','VerticalAlignment','middle', 'HorizontalAlignment','left', args{:});
v(6)=text( 0, 0, -1.1,'Dip 0','VerticalAlignment','middle', 'HorizontalAlignment','left', args{:});
axis([-1.3, 1.3 -1.2 1.2 -1.2 1.2])
hold off
title('Unitvectors coloured', 'fontsize', 14)
end
