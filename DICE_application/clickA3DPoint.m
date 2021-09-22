function h = clickA3DPoint(pointCloud)
%CLICKA3DPOINT
%   H = CLICKA3DPOINT(POINTCLOUD) shows a 3D point cloud and lets the user
%   select points by clicking on them. The selected point is highlighted 
%   and its index in the point cloud will is printed on the screen. 
%   POINTCLOUD should be a 3*N matrix, represending N 3D points. 
%   Handle to the figure is returned.
%
%   other functions required:
%       CALLBACKCLICK3DPOINT  mouse click callback function
%       ROWNORM returns norms of each row of a matrix
%       
%   To test this function ... 
%       pointCloud = rand(3,100)*100;
%       h = clickA3DPoint(pointCloud);
% 
%       now rotate or move the point cloud and try it again.
%       (on the figure View menu, turn the Camera Toolbar on, ...)
%
%   To turn off the callback ...
%       set(h, 'WindowButtonDownFcn',''); 
%
%   by Babak Taati
%   http://rcvlab.ece.queensu.ca/~taatib
%   Robotics and Computer Vision Laboratory (RCVLab)
%   Queen's University
%   May 4, 2005 
%   revised Oct 30, 2007
%   revised May 19, 2009

if nargin ~= 1
    error('Requires one input arguments.')
end

% if size(pointCloud, 1)~=3
%     error('Input point cloud must be a 3*N matrix.');
% end

% show the point cloud
h = gcf;
%pcshow(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:),pointCloud(:,4:6)/255);% !!!!!!!!!!!!!!!!!!
scatter3(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:), 3,pointCloud(4:6,:)'/255);
grid on
axis equal
% plot3(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:), 'c.'); 

PCview = fit3dplane_2(pointCloud(1:3,:)');%fit a 3D plane to the pointcloud to find the best view for point selection

if PCview(2) < 180
    azView = PCview(2)-180;
elseif PCview(2) > 180
    azView = 360-PCview(2)+180;
elseif PCview(2) == 180
    azView =0;
elseif PCview(2) == 0 || PCview(2) == 360
    azView = 180;
end
elView=90-PCview(1);
view(azView,elView)

cameratoolbar('Show'); % show the camera toolbar
hold on; % so we can highlight clicked points without clearing the figure

% set the callback, pass pointCloud to the callback function
set(h, 'WindowButtonDownFcn', {@callbackClickA3DPoint, pointCloud(1:3,:)}); 

function [exportPoint]=callbackClickA3DPoint(src, eventData, pointCloud)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009

point = get(gca, 'CurrentPoint'); % mouse click position
camPos = get(gca, 'CameraPosition'); % camera position
camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to

camDir = camPos - camTgt; % camera direction
camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
zAxis = camDir/norm(camDir);    
upAxis = camUpVect/norm(camUpVect); 
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % view rotation 

% the point cloud represented in the view frame
rotatedPointCloud = rot * pointCloud(1:3,:); 

% the clicked point represented in the view frame
rotatedPointFront = rot * point' ;

% find the nearest neighbour to the clicked point 
pointCloudIndex = dsearchn(rotatedPointCloud(1:2,:)', ... 
    rotatedPointFront(1:2));

h = findobj(gca,'Tag','pt'); % try to find the old point
selectedPoint = pointCloud(:, pointCloudIndex); 

if isempty(h) % if it's the first click (i.e. no previous point to delete)
    
    % highlight the selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'r.', 'MarkerSize', 20); 
    set(h,'Tag','pt'); % set its Tag property for later use   

else % if it is not the first click

    delete(h); % delete the previously selected point
    
    % highlight the newly selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'r.', 'MarkerSize', 20);  
    set(h,'Tag','pt');  % set its Tag property for later use

end
exportPoint = [selectedPoint(1,:),selectedPoint(2,:),selectedPoint(3,:)];
%fprintf('you clicked on point number %d\n', pointCloudIndex);

