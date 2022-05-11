function [dipdireigv, dipeigv]=Eigenvector(dipdir,dip,idx3)
% Eigenvector.m calculates mean vectos for a cluster with eigenvectors
%
%
% 12.10.2010 M. Markovaara-Koivisto, Aalto University School of Science and
% Technology, Finland

for n=1:length(dip);
    if dipdir(n)>180
        X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
        X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
        X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
    elseif dipdir(n)<=180
        X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
        X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
        X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
    end
end


%EIGENVALUES AND EIGENVECTORS
% The mean orientation of a cluster of orientation is the eigenvector
% associated with the highest eigenvalue

% Calculates products from direction cosines 
for n1=1:max(idx3)
xx(n1)=sum(X1b(find(idx3==n1)).^2);
xy(n1)=sum(X1b(find(idx3==n1)).*X2b(find(idx3==n1)));
xz(n1)=sum(X1b(find(idx3==n1)).*X3b(find(idx3==n1)));
yy(n1)=sum(X2b(find(idx3==n1)).^2);
yz(n1)=sum(X2b(find(idx3==n1)).*X3b(find(idx3==n1)));
zz(n1)=sum(X3b(find(idx3==n1)).^2);

% Orientation matrix
S(:,:,n1)=[xx(n1) xy(n1) xz(n1); xy(n1) yy(n1) yz(n1);xz(n1) yz(n1) zz(n1)];

% Orientation matrix normalized with number of obervation in the cluster
Sn=S./sum(idx3==n1);

% B eigenvectors, K eigenvalues
[Bi,Ki]=eig(Sn(:,:,n1));
B(:,:,n1)=Bi;
K(:,:,n1)=Ki;

% Eigenvector associated with the highest eigenvalue-> mean vector of the
% group of N vectors
Kmax(:,:,n1)=max(K(:,:,n1));
n2=find(Kmax(:,:,n1)==max(Kmax(:,:,n1)));
mean_vector(:,n1)= B(:,n2,n1);

% Eigenvector's x,y and z coordinates 
xeig(n1)=mean_vector(1,n1);
yeig(n1)=mean_vector(2,n1);
zeig(n1)=mean_vector(3,n1);

% Transforms the coordinates into dip/dipdir form
% Pole vector's orientations in quarters
if xeig(n1)>0 & yeig(n1)>0
dipeigv(n1)=90-abs((asin(zeig(n1))))*(180/pi);
dipdireigv(n1)=abs((atan(xeig(n1)/yeig(n1))))*180/pi;
elseif xeig(n1)>0 & yeig(n1)<0
dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
dipdireigv(n1)=180-abs((atan(xeig(n1)/yeig(n1))))*180/pi;
elseif xeig(n1)<0 & yeig(n1)<0
dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
dipdireigv(n1)=180+abs((atan(xeig(n1)/yeig(n1))))*180/pi;
else
dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
dipdireigv(n1)=360-abs((atan(xeig(n1)/yeig(n1))))*180/pi;
end
end

% Prints number of observations in the cluster and dip/dip direction into
% the RESULTS structure
for n1=1:max(idx3)
   RESULTS.(genvarname(['Set', num2str(n1)])).ObservationNumberInCluster=sum(idx3==n1);
   RESULTS.(genvarname(['Set', num2str(n1)])).Dip=dipeigv(n1);     
   RESULTS.(genvarname(['Set', num2str(n1)])).Dipdir=dipdireigv(n1);
end
end
