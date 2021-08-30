function P32calculators(pn,fn, PC, Fracdata)
%Function that calculate the P32 using the fracture geometry dataset and
%the pointcloud;
%This function calculate the P32 for each point of the point cloud

tic
% Do you want to consider the sets ? yes q=1, no q=0
q=0;

%% 1) LOAD & READ XLSX CIRCULAR FRACTURE DATA
%  Select CSV file containing all plane
if ismac %For OSX operative system
    [fn, pn] = uigetfile({'*.xlsx', 'Select a XLSX fracture geometry file'},'Select a fracture geometry file');% <- MODIFY the PATH
elseif ispc % For WINDOWS operative system
    [fn, pn] = uigetfile({'*.xlsx', 'Select a XLSX fracture geometry file'},'Select a fracture geometry file');% <- MODIFY the PATH
end
Fracdata=readtable(fullfile(pn, fn));
nplane=numel(Fracdata.Dip); %number of rows (= number of plane)
%% 2) LOAD & READ XLSX SET BELONGING DATA
if q==1
    [fneset, pnset] = uigetfile({'*.xlsx', 'Select a XLSX fracture set file'},'mytitle',...
        pn);
    Fracset=readtable(fullfile(pnset, fnset));
    nplane=numel(Fracset.Set); %number of rows (= number of plane)
end
%% LOAD & READ the Pointcloud file
usePC=1; %if UsePC=1 the code use the pointcloud, if = 0 it doesn't use the PC
[fnPC, pnPC] = uigetfile({'*.txt', 'Select a pointcloud in txt'},'Select a pointcloud',pn);
PC=importdata(fullfile(pnPC, fnPC));
%% 4) DEFINE VARIABLES FOR CALCULATION
dipdir =  Fracdata.DipDirection(:);%Dip direction
dip =  Fracdata.Dip(:);%Dip angle
xyz = [Fracdata.Xcenter(:), Fracdata.Ycenter(:), Fracdata.Zcenter(:)];
Nxyz = [Fracdata.Nx(:),Fracdata.Ny(:),Fracdata.Nz(:)];
radius = Fracdata.Radius(:);

if q==1
    Set = Fracset.Set;
    Set(isnan(Set)) = max(Set)+1;
    Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};% define color for set from 1 to 7
end
if usePC==0
    SVxyz=[1.295830,1.103351,1.608449];% define scan volume (sphere) center
    PC=SVxyz;
end

% define scan volume (sphere) radius


%%P32 calculation
%1) calculate intersection between the scan volume and the fracture
%'infinite plane'
sumarea=0;
if q==1
    sumareaset=zeroz(1,max(Set));
end
SVrrange = 0.1:0.1:1; %range of radius of the scan volume

SVresults=zeros(length(SVrrange),2);
P32=zeros(length(PC),1);
mkdir (pn, 'exported_pointclouds_parfor')
disp('DATASET_STATISTIC')
disp(['Pointcloud_= ', num2str(length(PC(:,1))),' (points)'])
disp(['Fracture_numbers_= ', num2str(nplane)])
%% CALCULATE P32
for k = 1 : length(SVrrange)
    SVr=SVrrange(k);
    %rrange(k);
    clear P32;
    clear P32set;
    clear tic
    disp(['START_calculation_P32_for_scan radius_= ',num2str(SVr,'%.1f'),' (m)'])
    tic
    parfor j=1:length(PC(:,1))
        
        % % % %         fid = fopen(fullfile([pn,'exported_pointclouds_parfor'],['P32_r_',num2str(SVr,'%.1f'),'.txt']),'wt');
        SVxyz= PC(j,1:3);
        SV=[SVxyz,SVr];% scan volume express as required yby the function intersectPlaneSphere
        sumarea=0;
        
        for i=1 : nplane
            %PLANEtemp=zeros(1,6);
            scanCIRCtemp =zeros(1,4);
            %IntArea=0;
            
            PLANEtemp = [xyz(i,:), Nxyz(i,:)];
            
            %dPlCi = bsxfun(@minus, xyz(i), SVxyz); %distance between the centers of the plane and the sphere
            %dPlCi = sqrt(d(1)^2 + d(2)^2 + d(3)^2);
            dOrSp= (-dot (Nxyz(i,:), SVxyz)); %distance from origin-sphere centertoward normal vector of the plane
            dOrPl= (-dot (Nxyz(i,:), xyz(i,:))); %distance from origin-plane centertoward normal vector of the plane
            aSpPl=abs(dOrSp - dOrPl); %distance from sphere-plane centers toward normal vector of the plane
            
            
            
            if aSpPl < SVr % && dPlCi<(SVr+radius(i))
                scanCIRCtemp = intPlaneSphere(PLANEtemp, SV);
                IntArea = intCirc2Circ(scanCIRCtemp(1:3),scanCIRCtemp(4), xyz(i,:), radius(i));
            else
                IntArea = 0;
            end
            
            
            
            sumarea = sumarea + IntArea;
        end
        
        P32(j,1) = sumarea / ((4/3)*pi*SVr^3); %total P32 (all sets together)
        % % %         if q==1;
        % % %             P32set = sumareaset / ((4/3)*pi*SVr^2);%P32 specific for each set)
        % % %         end
        % % %
        % % %         %P32(isnan(P32))=0;
        % % %         while fid == -1
        % % %             pause(0.05);
        % % %             fid = fopen(fullfile([pn,'exported_pointclouds_parfor'],['P32_r_',num2str(SVr,'%.1f'),'.txt']),'wt');
        % % %         end
        % % %         fprintf(fid, '%.3f %.3f %.3f %.3f \n', [SVxyz,P32]);
        % % %         fclose(fid);
    end
    
    
    %SVresults(k,1)=rrange(k);
    %SVresults(k,2)=P32;
    %fclose(fid);
    disp(['->Exporting_results_at_= ',num2str(toc/60,'%.1f'), '(min)'])
    clear fid
    fid = fopen(fullfile([pn,'exported_pointclouds_parfor'],['P32_r_',num2str(SVr,'%.1f'),'.txt']),'wt');
    fprintf(fid, 'X Y Z R G B P32 \n');
    fprintf(fid, '%.6f %.6f %.6f %3.0f %3.0f %3.0f %.3f \n', [PC(:,1),PC(:,2),PC(:,3),PC(:,4),PC(:,5),PC(:,6),P32(:)]');
    fclose(fid);
    % Consider to use dlmwrite
    % dlmwrite(fullfile([pn,'exported_pointclouds_parfor'],['P32_r_',num2str(SVr,'%.1f'),'.txt'])',...
    % [PC(:,:),P32(:)],'precision','%.3f','delimiter',' ');
    disp('End of point cloud exporting process')
end
end

