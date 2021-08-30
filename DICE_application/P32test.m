close all
clear variable
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

SVr=1.0;% define scan volume (sphere) radius


%%P32 calculation
%1) calculate intersection between the scan volume and the fracture
%'infinite plane'
sumarea=0;
if q==1
    sumareaset=zeroz(1,max(Set));
end
rrange = 1;%0.1:0.1:10; %range of radius of the scan volume

SVresults=zeros(length(rrange),2);

mkdir (pn, 'exported_pointclouds')
%% CALCULATE P32
for k = 1 : length(rrange)
    fid = fopen(fullfile([pn,'exported_pointclouds'],['P32_r_',num2str(rrange(k),'%.1f'),'.txt']),'wt');
    SVr=rrange(k);
    clear P32;
    clear P32set;
    for j=1:length(PC(:,1))
        clear SVxyz
        
        SVxyz=PC(j,1:3);
        SV=[SVxyz,SVr];% scan volume express as required yby the function intersectPlaneSphere
        
        
        sumarea=0;
        if q==1
            sumareaset=zeroz(1,max(Set));
        end
        
        for i=1 : nplane
            clear PLANEtemp CIRCtemp IntArea;
            PLANEtemp = [xyz(i,:), Nxyz(i,:)];
            
            dPlCi = bsxfun(@minus, xyz(i), SVxyz); %distance between the centers of the plane and the sphere
            dPlCi  = sqrt(dPlCi(1)^2 + dPlCi(2)^2 + dPlCi(3)^2);
            
            dOrSp= (-dot (Nxyz(i,:), SVxyz)); %distance from origin-sphere centertoward normal vector of the plane
            dOrPl= (-dot (Nxyz(i,:), xyz(i,:))); %distance from origin-plane centertoward normal vector of the plane
            aOrPl=abs(dOrSp - dOrPl); %distance from sphere-plane centers toward normal vector of the plane
            
            
            
            if aOrPl < SVr % && dPlCi<(SVr+radius(i))
                CIRCtemp = intPlaneSphere(PLANEtemp, SV);
                IntArea = intCirc2Circ(CIRCtemp(1:3),CIRCtemp(4), xyz(i,:), radius(i));
            else
                IntArea = 0;
            end
            
            
            
            sumarea = sumarea + IntArea;
            
            if q==1
                sumareaset(1,Set(i)) = sumareaset(1,Set(i)) + IntArea;
            end
            
            
            
            
        end
        
        P32 = sumarea / ((4/3)*pi*SVr^2); %total P32 (all sets together)
        if q==1
            P32set = sumareaset / ((4/3)*pi*SVr^2);%P32 specific for each set)
        end
        
        %P32(isnan(P32))=0;
        fprintf(fid, '%.3f %.3f %.3f %.3f \n', [SVxyz,P32]);  
    end
    SVresults(k,1)=rrange(k);
    SVresults(k,2)=P32;
    fclose(fid);
    
end
toc