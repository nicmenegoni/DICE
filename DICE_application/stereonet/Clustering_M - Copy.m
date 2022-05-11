
% Clustering_M.m clusters the data around pre-defined centroids with kmeans-
% method.
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland

m=nplane;
% Clustering is carried out on the discontinuities polevectors, which are 
% calculated into N-matrice.

% Calculation of the surface's normal
% Surface is defined with two perpendicular vectors V and O
V=zeros(1,3,m);
for r=1:m
% Unit vector V, which origin is at the origo, in the orientation of the 
% dip/dip direction
V(:,:,r)=[cosd(90-dipdir(r))*cosd(dip(r)) sind(90-dipdir(r))*cosd(dip(r)) -sind(dip(r))];
% Unit vector O, which origin is at the origo, which is perpendicular to 
% the dip direction and is horizontal
O(:,:,r)=[-sind(90-dipdir(r)) cosd(90-dipdir(r)) 0];
end
% Normal vector N of the surface defined with the vectors V and O
N=zeros(1,3,m);
for k=1:m 
    N(:,:,k)=cross(V(:,:,k),O(:,:,k));
end
% Pole vector is counter vector of the normal vector, and it points
% downwards 
N=-N;
antiN=-N;



% Asks to klick on the approximate centroids of the clusters nclu times
disp('Left click with mouse on the approximate centroids of the clusters.')
[xclu yclu]=ginput(nclu);

% Middlepoints of the clusters
% Clusters' coordinates are transformed into Dip/Dipdirection form
cludipdir=zeros(1,nclu);
cludip=zeros(1,nclu);

for s=1:nclu
if xclu(s)>0 & yclu(s)>0
    clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
    cludipdir(s)=270-clutheta(s)*180/pi;
    cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
    cludip(s)=abs(cludipr(s)*180/pi);
elseif xclu(s)>0 & yclu(s)<0
    clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
    cludipdir(s)=270+clutheta(s)*180/pi;
    cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
    cludip(s)=abs(cludipr(s)*180/pi);
elseif xclu(s)<0 & yclu(s)<0
    clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
    cludipdir(s)=90-clutheta(s)*180/pi;
    cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
    cludip(s)=abs(cludipr(s)*180/pi);
else
    clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
    cludipdir(s)=90+clutheta(s)*180/pi;
    cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
    cludip(s)=abs(cludipr(s)*180/pi);
end
end

% Dip/Dip direction forms is transformed into polevectors, which are used in the clustering
% Calculation of the normal vectors N
% Defining the surface with two perpendicular vectors V and O
V=zeros(1,3,nclu);
for r=1:nclu
% Unit vector V, which origin is at the origo, in the orientation of the 
% dip/dip direction
V(:,:,r)=[cosd(90-cludipdir(r))*cosd(cludip(r)) sind(90-cludipdir(r))*cosd(cludip(r)) -sind(cludip(r))];
% Unit vector O, which origin is at the origo, which is perpendicular to 
% the dip diprecion and is horizontal
O(:,:,r)=[-sind(90-cludipdir(r)) cosd(90-cludipdir(r)) 0];
end
% Normal vector N of the surface defined with the vectors V and O
Nclu=zeros(1,3,s);
for k=1:s 
    Nclu(:,:,k)=cross(V(:,:,k),O(:,:,k));
end
%Pole vector is counter vector of the normal vector, and it points
%downwards 
Nclu=-Nclu;



% Iteration of the contents of the clusters: Initial centroids of the
% clustres are given and the new ones are calculated, observations are moved
% to the clusters and sum of distanced of centroids and the observations 
% within the cluster are calculated
% First iteration
previous_clustercosum=[];
% Observation are divided into the closest clusters
arcdistance=zeros(s,m);
for s=1:nclu
    for k=1:m
        arcdistance(s,k)=(acos(dot(N(:,:,k),Nclu(:,:,s))))^2; % distance between the observation and the centroids, nclu defines rows and number of observations columns
        antiarcdistance(s,k)=(acos(dot(antiN(:,:,k),Nclu(:,:,s))))^2;
    end
end
for s=1:nclu
    for k=1:m
        rownumber(k)=find(arcdistance(:,k)==min(arcdistance(:,k))); % Searches the rownumber of the centroid, which has the smallest distance to the observation point.
        idx3(k)=rownumber(k); % Observation gets the rownumber as its cluster number.
    end
end
idx3first=idx3;
% Calculates the sum of arcdistances within the cluster
for s=1:nclu
        clustersum(s)=sum(arcdistance(find(idx3==s)));
end
% Sums up all the cluster's sums together
cosum=sum(clustersum);

% Iterates the cluster cosum as long as it does not get any smaller

% Second and later iterations    
while (isempty(previous_clustercosum) | cosum < previous_clustercosum)
    % Calculates new centroids of the clusters
    % Calculates sumvector of the clusters' polevectors
    previous_clustercosum=cosum;
    for s=1:nclu
            sum_cluster(:,s)=[sum(N(1,1,find(idx3==s))) sum(N(1,2,find(idx3==s))) sum(N(1,3,find(idx3==s)))]; %Finds the rows with the same idc3 and sums them up
            % Length of the polevectors' sum
            L(s)=sqrt(sum(N(1,1,find(idx3==s)))^2+sum(N(1,2,find(idx3==s)))^2+sum(N(1,3,find(idx3==s)))^2);
            % Sumvector is normalised into a unitvector
            new_Nclu(:,s)=sum_cluster(:,s)/L(s);
    end

    % Points are re-divided according to the closest new cluster centroids
    arcdistance=zeros(s,m);
    for s=1:nclu
        for k=1:m
           arcdistance(s,k)=(acos(dot(N(:,:,k),new_Nclu(:,s))))^2; % distance between the observation and the centroids, nclu defines rows and number of observations columns
           
        end
    end
    for s=1:nclu
        for k=1:m
        rownumber(k)=find(arcdistance(:,k)==min(arcdistance(:,k))); % Searches the rownumber of the centroid, which has the smallest distance to the observation point.
        end
    end
    
    idx3=rownumber; % Observation gets the rownumber as its cluster number.
    idx3second=idx3;
    % Calculates the sum of arcdistances within the cluster
    for s=1:nclu
        clustersum(s)=sum(arcdistance(find(idx3==s)));
    end
    % Sums up all the clustre's sums together
    cosum=sum(clustersum);
    disp('### iteration ###')
end
