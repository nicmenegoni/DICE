clear all
close all
% function [c_Poly,c_Poi] = f_LectDxf(nomArch)
% FUNCTION dxf = read_dxf(filename)
%
% Author:  Steven Michael (smichael@ll.mit.edu)
%
% Date:    3/10/2005
%
% Description
%
%   The following compiled MATLAB file reads in an
%   ascii DXF file and loads the information into
%   the "dxf" variable
%
% Inputs:
%
%   Filename   :    The filename of the ASCII DXF File
%
% Outputs:
%
%   dxf        :    A 3-D variable of size (NX3X3).  The first
%                   index N is the number of facets.  The second
%                   index references the 3 vertices of each
%                   facet, and the third index is the
%                   (x,y,z) location of the vertex.
%
% Note:
%
%   This may not work with all DXF files.  The DXF file format
%   is complicated and continuously evolving.  This has worked
%   with every file I happened to use it on, but that will
%   not be true for everyone.
%
%   Also, the function does not load any color or texture
%   information.
%%% Read entities information of dxf file
%author = SebaTM
%Jul 27 2009
%Based in dxf2coord 1.1 matrix of lukas wischounig, but is not dependent of the Code
%Group relative position. That is better way to read dxf format. Don't fail if the
%polyline has arcs (budges), but yet don't read them. Don't read arcs as circles. Read
%properties (see case 'LINE' by examples of modifications). Group Codes and Associated
%Values is read in accurately way (accept associated values with space).
%
%Use cell2mat(cell(:,1)) to acquire geometry data in matrix,
%by example cell2mat(c_Cir(:,1))
clear all
close all
[filename, pathname]=uigetfile({'*.dxf', 'Select a DXF file'}, 'Select a DXF file',...
    'F:\Menegoni\Antola\Outcrop_models');
    %'F:\Menegoni\Antola\Outcrop_models');% <- MODIFY the PATH
% select path in which DXF of Circular Plane will be saved
disp('########### START OF FITTING PROCCES ##########')
tic
%% Read file
fId = fopen(fullfile(pathname,filename));
c_ValAsoc = textscan(fId,'%d%s','Delimiter','\n');
fclose(fId);
% Code Group Matrix
m_GrCode = c_ValAsoc{1};
% Associated value String Cell
c_ValAsoc = c_ValAsoc{2};
%[m_GrCode,c_ValAsoc] = c_ValAsoc{:};

%% Entities
m_PosCero = find(m_GrCode==0);
%Is searched by (0,SECTION),(2,ENTITIES)
indInSecEnt = strmatch('ENTITIES',c_ValAsoc(m_PosCero(1:end-1)+1),'exact');
%(0,ENDSEC)
m_indFinSecEnt = strmatch('ENDSEC',c_ValAsoc(m_PosCero(indInSecEnt:end)),'exact');
% Entities Position
m_PosCero = m_PosCero(indInSecEnt:indInSecEnt-1+m_indFinSecEnt(1));
% Variable initiation
%accelerate?
     c_Line = cell(sum(strcmp('LINE',c_ValAsoc(m_PosCero))),2);
     c_Poly = cell(sum(strcmp('POLYLINE',c_ValAsoc(m_PosCero))),2);
     c_Cir = cell(sum(strcmp('CIRCLE',c_ValAsoc(m_PosCero))),2);
     c_Arc = cell(sum(strcmp('ARC',c_ValAsoc(m_PosCero))),2);
     c_Poi = cell(sum(strcmp('POINT',c_ValAsoc(m_PosCero))),2);
c_Line = cell(1,2);
c_Poly = cell(1,2);
c_Cir = cell(1,2);
c_Arc = cell(1,2);
c_Poi = cell(1,2);
%
iLine = 1;
iPoly = 1;
iCir = 1;
iArc = 1;
iPoi = 1;
% Loop on the Entities
for iEnt = 1:length(m_PosCero)-2
    m_GrCodeEnt = m_GrCode(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    c_ValAsocEnt = c_ValAsoc(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    nomEnt = c_ValAsocEnt{1};  %c_ValAsocEnt{m_PosCero(iEnt+1)}
    %In the entitie's name is assumed uppercase
    switch nomEnt
        case 'VERTEX'
            % (X,Y,Z) Position
            c_Poi{iPoi,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                str2double(f_ValGrCode(30,m_GrCodeEnt,c_ValAsocEnt))];
            % Layer
            c_Poi(iPoi,2) = f_ValGrCode(8,m_GrCodeEnt,c_ValAsocEnt);
            % Add properties
            %
            iPoi = iPoi+1;
            %case Add Entities
        
            % (X,Y,Z) vertexs
            
   
    end
end
%% create XYZNCloud matrix with X, Y and Z coordinate and Number of pointcloud
for i=1:numel(c_Poi)/2
    
    clear temp_name
    temp_name=char(c_Poi(i,2));
    Ncloud_c_Poi (i,1) = str2num(temp_name(10:end));
    
end
XYZNCloud=cell2mat(c_Poi(:,1));
XYZNCloud(:,4)=Ncloud_c_Poi(:,1);
%end
%%
figure(1)
plungeRad=zeros(1,max(XYZNCloud(:,4)));
for i=1 : max(XYZNCloud(:,4))%max(XYZNCloud(:,4)=number of planes
%for n= 1, 2
%[x, y, z,] = textread('/Users/niccolomenegoni/Desktop/DATI Windows/università/Dottorato Pavia/DATI_PhD/Antola/St_280116/Stazione 3bis_new/lines/lines.txt', ...
%'%f %f %f', 1)
 %St5
 clearvars xyztemp
 xyztemp=XYZNCloud(XYZNCloud(:,4)==i,:);
 x1(i)=xyztemp(1,1);
 y1(i)=xyztemp(1,2);
 z1(i)=xyztemp(1,3);
 
 x2(i)=xyztemp(2,1);
 y2(i)=xyztemp(2,2);
 z2(i)=xyztemp(2,3);
 
%  
 L(i)= sqrt((x2(i)-x1(i))^2+(y2(i)-y1(i))^2+(z2(i) -z1(i))^2);
 cosAlpha(i)=(x2(i)-x1(i))/L(i);
 cosBeta(i)=(y2(i)-y1(i))/L(i);
 cosGamma(i)=(z1(i)-z2(i))/L(i);
 plungeRad(i) = asin(- cosGamma(i));
 azimuthRad(i) = atan((x2(i)-x1(i))/(y2(i)-y1(i))); 
 
 plungeDeg(i)=radtodeg(plungeRad(i));
 azimuthDeg(i)=radtodeg(azimuthRad(i));
 
 
 
 %disp('########### AZIMUTH,PLUNGE CORRECTED ##############')
 if cosBeta(i)<0
     azimuthDeg_c(i)= 180 + azimuthDeg(i);
 else
     if cosAlpha(i)>0
         azimuthDeg_c(i)= azimuthDeg(i);
     else cosAlpha(i)<0
         azimuthDeg_c(i)= 360+azimuthDeg(i);
     end
 end
 
 %disp('###########  Corrected ##############')
 if plungeDeg(i)<0
     plungeDeg_c(i)=-plungeDeg(i);
     azimuthDeg_c(i)=azimuthDeg_c(i);
 else
     plungeDeg_c(i)=plungeDeg(i);
     if azimuthDeg_c(i) > 180 || azimuthDeg_c(i)==0
         azimuthDeg_c(i)=azimuthDeg_c(i)-180;
     else
         azimuthDeg_c(i)=azimuthDeg_c(i)+180;
     end
 end
  azimuthRad_c(i)=degtorad(azimuthDeg_c(i));
  plot3([x1(i),x2(i)],[y1(i),y2(i)],[z1(i),z2(i)])
  hold on
  axis equal
  grid on
end
figure(1)
title(['N = ',num2str(length(L))])
xlabel('East (m)')
ylabel('North (m)')
zlabel('Elevation (m)')
savefig(fullfile(pathname,'3Dplot'))
 figure(2)
 title(['N = ',num2str(length(L))])
 polarhistogram(azimuthRad_c,36);
 set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise')
 figure(2)
savefig(fullfile(pathname,'RosePlot'))
%Export values
Values=zeros(length(L),9);
Values(:,1)=round(azimuthDeg_c(:));
Values(:,2)=round(plungeDeg_c(:));
Values(:,3)=L(:);
Values(:,4)=x1(:);
Values(:,5)=y1(:);
Values(:,6)=z1(:);
Values(:,7)=x2(:);
Values(:,8)=y2(:);
Values(:,9)=z2(:);
varname(1,1)={'Azimuth_deg'};
varname(1,2)={'Plunge_deg'};
varname(1,3)={'Length_m'};
varname(1,4)={'x1'};
varname(1,5)={'y1'};
varname(1,6)={'z1'};
varname(1,7)={'x2'};
varname(1,8)={'y2'};
varname(1,9)={'z2'};
Table=array2table(Values);
Table.Properties.VariableNames=varname;
tablefilenameXLSX='AzimuthPlungeResults.xlsx';
writetable(Table,fullfile(pathname,tablefilenameXLSX),'WriteRowNames',true);
