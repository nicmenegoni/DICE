function plane3d = fit3dplane(PC)
% This function is developed following the work of Jones et al. (2015)
% and Seers and Hodgetts (2016)
% 
% INPUT
% 'PC' (point cloud) is a vector or a matrix N x 3 (N=number of points)
% that contain X Y Z coordinate of the points of the cloud

centroid = mean(PC(:,1:3)); % calculate centroid of point cloud

%[U, S, W] = svd(PC);

C = cov(PC); % calculate the covariance matrix 'C'

[V, D] = eig(cov(PC)); % calculate eigenvector 'V' and eigenvalue 'D' of the
% covariance matrix 'C'


M = log(D(3,3)/D(1,1));% calculate vertex coplanarity 'M'
K = log(D(3,3)/D(2,2))/log(D(2,2)/D(1,1));% calculate vertex collinearity 
%'K'


N = V(:,1)';

% Normals have to be oriented upward!!! (Normals have to be horiented in same
% direction)

if N(1,3)<0 %if normal is oriented downward
    N(1,1)=-N(1,1);
    N(1,2)=-N(1,2);
    N(1,3)=-N(1,3);    
end
cosAlpha = N(1,1)/norm(N);
cosBeta = N(1,2)/norm(N);
cosGamma = N(1,3)/norm(N);

dip = 90 + rad2deg(asin(-cosGamma));

if cosAlpha > 0 && cosBeta > 0
dipdir = rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha > 0 && cosBeta < 0
dipdir = 180 + rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha < 0 && cosBeta < 0
dipdir = 180 + rad2deg(atan(cosAlpha/cosBeta));
elseif cosAlpha < 0 && cosBeta > 0
dipdir = 360 + rad2deg(atan(cosAlpha/cosBeta));
end

%calculate radius of plane = max distance point-centroid of fitted plane
for i=1:numel(PC(:,1))
pdist1 (i) =sqrt((centroid(1,1)-PC(i,1))^2 + (centroid(1,2)-PC(i,2))^2 + (centroid(1,3)-PC(i,3))^2);
end

[pM1, pI1] = max(pdist1);

for i=1:numel(PC(:,1))
    if i==pI1
        pdist2(i)=0;
    else
        pdist2 (i) =sqrt((PC(pI1,1)-PC(i,1))^2 + (PC(pI1,2)-PC(i,2))^2 + (PC(pI1,3)-PC(i,3))^2);
    end
end
[pM2, pI2] = max(pdist2);

% centroid=mean([PC(pI1,:); PC(pI2,:)]);
radius =(sqrt((PC(pI1,1)-PC(pI2,1))^2 + (PC(pI1,2)-PC(pI2,2))^2 + (PC(pI1,3)-PC(pI2,3))^2))/2;

%calculation of mean error and standard deviation of fitted plane
d = -sum(bsxfun(@times, N, bsxfun(@minus, centroid, PC)), 2);
mean_err=mean(abs(d));
stdev_err=std(abs(d));


plane3d = [dip, dipdir, radius, centroid, N, M, K, mean_err, stdev_err];


% REFERENCES:
%
% - Jones, R. R., Pearce, M. A., Jacquemyn, C., & Watson, F. E. (2016). 
%   Robust best-fit planes from geospatial data. Geosphere, 12(1), 196-202.
%
% - Seers, T. D., & Hodgetts, D. (2016). Probabilistic constraints on
%   structural lineament best fit plane precision obtained through
%   numerical analysis. Journal of Structural Geology, 82, 37-47.