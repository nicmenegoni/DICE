close all
clear all


[filename, pathname] = uigetfile({'*.csv', 'Select a CSV file'},'mytitle',...
          'F:\DATI\D_data\dottorato\DATI\Antola\Outcrop_models');
if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
end
%reading text file in which coloumn are: 1,2,3=x,y,z and 4,5,6=Nx,Ny,Nz
% filename = '/Users/niccolomenegoni/Desktop/220_70.txt'
[CSV,delimiterOut]=importdata(fullfile(pathname, filename),';',1);
n=numel(CSV.data)/16; %number of rows
disp ([filename])
A=zeros(n ,6);

for i= 1 : n
    for j = 1 : 6
        
        A(i,j)= CSV.data(i,j+1);
      
    end
end

dipdir =  CSV.data(:,13);%Dip direction
dip =  CSV.data(:,14);%Dip angle
nplane=n; % number of plane to plot in the stereogram
[xp, yp]=Stereogram_colourless_Schmidt(dipdir, dip, nplane); %%%% FINCO A QUA CORRETTO




% Draws the oriented data on a Schmidt's equal are net


% Clustering of the discontintuities
% Asks the number of joint sets in the data
nclu = input('Number of joint sets (+1 for randoms) in your data (max 7):');
Clustering_M

% Draws the automated clustering results (kmeans-method) on a Schmidt's equal area net
 Stereogram(dipdir, dip, nplane, idx3, xp, yp)

% Asks is the user is satisfied with the result of the clustering
% If not, the user may change the clustering in an interactive window
Ready = [];
button=1;
answer = input('Do you want to change the results of the automatic clustering (y/n):', 's');
 while (answer == 'y' & isempty(Ready) & button==1)
 disp('Left click with mouse on the symbols as many times as it takes to change the cluster into the wanted one.')
 disp('When you are finnished, right click with mouse on the figure')
 [xcho, ycho, button]=ginput(1);

% Calculates the minimun distance between the indicated point and the data
% points and searches the closes data point, gives index to the row and 
% changes the cluster group into the next one or to the first one.
% Refreshes the figure.
 distance=sqrt((xp-xcho).^2+(yp-ycho).^2);
 cho=find(distance==min(distance))
 % Observation is indicated to the next cluster
    if (idx3(cho)<max(idx3) & button==1)
     idx3(cho)=idx3(cho)+1;
    elseif button==1
        idx3(cho)=1;
    end
 Stereogram
 end

 % Draws colour coded unit vectors on a stereographic sphere, in which dip 0
 % degrees is at the bottom of the sphere and 90 degrees at the top.
 Unitvectors_coloured(dipdir, dip, nplane, idx3)
 
 % Calculates the mean attitudes of the clusters
[dipdireigv, dipeigv]=Eigenvector(dipdir, dip, idx3);
 
 % Draws a stereographic projection with the mean attitudes of the clusters
StereogramEig(dipdir, dip, nplane, idx3)

    thetaeig=(90-dipdireigv);
    reig=sqrt(2)*sind((90-dipeigv-90)/2);
 for n=1:max(idx3)
    xeigp(n) = reig(n)*cosd(thetaeig(n));
    yeigp(n) = reig(n)*sind(thetaeig(n));
 end
 plot(xeigp, yeigp, 'kO')
