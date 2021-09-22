function [CWplane, maxCW_dist, minCW_dist] = selectCW(PC,CWr)
CWpoints=zeros(1,3);
Ready1 = [];button1=1;answer1 = 'y';
while (answer1 == 'y' & isempty(Ready1) & button1==1)
    clearvars CWpoints
    Ready2 = [];
    button2=1;
    answer2 = 'y';
    while (answer2 == 'y' & isempty(Ready2) & button2==1)
       %h = clickA3DPoint(PC(:,1:3)');
        h = clickA3DPoint(PC');
        answer2 = input('Do you want to change points? (y/n):', 's');
    end
    %[preCWxyz]=callbackClickA3DPoint(h, 'WindowButtonDownFcn', PC(:,1:3)');
    [preCWxyz]=callbackClickA3DPoint(h, 'WindowButtonDownFcn', PC');
    count=1;
    for i =1 : length(PC)
        d = sqrt((preCWxyz(1)-PC(i,1))^2 + (preCWxyz(2)-PC(i,2))^2 +...
            (preCWxyz(3)-PC(i,3))^2);% distance of point from the center of the sphere
        
        
        if d < CWr
            CWpoints(count,1)=PC(i,1);
            CWpoints(count,2)=PC(i,2);
            CWpoints(count,3)=PC(i,3);
            count=count+1;
        end
    end
    
    figure
    pcshow(PC(:,1:3),PC(:,4:6)/255);
    hold on
    pcshow(CWpoints,'y');
    CWplane = fit3dplane_2(CWpoints);
    
    theta=0:0.1:2*pi; %Setting the spacing in which coordinate of the discontinuity disc are calculated
    v=null(CWplane(7:9));% calculate vectors needed to plot the discs
    points=repmat(CWplane(4:6)',1,size(theta,2))+CWr*(v(:,1)*cos(theta)+v(:,2)*sin(theta));%calculate points coordinate of the discs
    hold on;grid on;box on
    fill3(points(1,:),points(2,:),points(3,:),'r', 'FaceAlpha', 0.5);%Plot 3D disc for each plane
    xlabel('x-axis (East)');ylabel('y-axis (Nord)');zlabel('z-axis (Elev)')
    
    answer1 = input('Do you want to discard this CW and select another:', 's');
end

for i = 1 : length(CWpoints)
    dn (i) = plane2point(CWplane(7:9),CWplane(4:6), CWpoints(i,:));
end
maxCW_dist = max(dn);
minCW_dist = min(dn);


end