function [L,l,m_mauldon,n_mauldon,countCWDisc_inter] = intCWdisc(CWxyz, CWNxyz, CWr, xyz, Nxyz, radius, Set, f1)
%% INFORMATION
% This function calculates the intersection between the defined circular
% window and the dataset of fractures/discontinuities and calcuates the
% P21 (areal intesity) and the fracture/discontinuity intersection traces
% statistic, as the Mean Trace Length (MTL) and the Trace length standard
% deviation.
% Moreover, the funtion use the correction method proposed by Mauldon
% (2001) to obtain the 'unbiased' value of Intesity (I), Density (d), Mean
% Trace length (mu)
% Input parameters:
%   CWxyz = center of the Circular Window;
%   CWNxyz = normal vector of the Circular Window;
%   CWr = radius of the Circular Window;

%% Plot stuff
%fit disc of circular window
nplane = length(radius);

%% Circular window - discontinuities intersections
% two arrays are initialized for the 'm' and 'n' Mauldon
% correction values
m_mauldon=zeros(nplane,1);% m Mauldon value
n_mauldon=zeros(nplane,1);% n Mauldon value

% initializig a matrix for counting the intersection of each discontinuity and set
countCWDisc_inter=zeros(nplane,max(Set));%
figure(f1)
hold on
grid on
% Preallocating memory
lineAB=zeros(nplane,6);
L=zeros(nplane,1);
l=zeros(nplane,1);
Point_intersCW_disc=zeros(nplane,6);

for i=1:nplane
    %%Compute intersection between two plane
    tol = 1e-14;% setting the angle cutoff for which pseudo parallel planes (parallel at the CW plane) are not consider.
    if abs(cross(Nxyz(i,:), CWNxyz, 2)) > tol
        [Pline, Vline] = intersect_planes(CWxyz, CWNxyz, xyz(i,:), Nxyz(i,:));
        lineAB(i, 1:6)=[Pline,Vline];%%Infinite line representing intersection between two infinte planes of discontinuity and CW
        
        point1AB(i,:)= xyz(i,:);% Center of discontinuity
        point2AB= CWxyz;% center of Circular Window
        
        %Distance between disconintuity center and previous intersection line
        dpointL_PA(i) = distPointLine(xyz(i,:),Pline,Vline);
        
        %Distance between circular window center and previous intersection line
        dpointL_PB(i) = distPointLine(CWxyz,Pline,Vline);
        
        if abs(dpointL_PA(i))<radius(i) && abs(dpointL_PB(i))<CWr %if distance is lower than radius of 2 discs (discont. and CW) intersection is considered
            clear d_point1A d_point1B d_point2A d_point2B d_point3A d_point3B d_point4A d_point4B
            selectAB=zeros(4,4);
            
            %Define Sphere (fitted for the 2 considered discs) for Sphere-Line Intersection
            SPHEREA = [xyz(i,1), xyz(i,2), xyz(i,3),  radius(i)];%definisco la sfera per la discontinuità
            SPHEREB = [CWxyz(1,1), CWxyz(1,2), CWxyz(1,3),  CWr];%definisco la sfera per la CW
            
            %Use Spehere-Line intersection function of geom3d package
            %(INRA) to find points of intersection
            PTSA= intersectLineSphere(lineAB(i,:), SPHEREA);% for line and discontinuity sphere
            P_interA(i,1:6)=[PTSA(1,1:3), PTSA(2,1:3)];
            PTSB = intersectLineSphere(lineAB(i,:), SPHEREB);% for line and CW sphere
            P_interB(i,1:6)=[PTSB(1,1:3), PTSB(2,1:3)];
            
            
            
            % Calculation of l, the cutted trace length, trace visble
            % inside the CW
            
            %calculation distance between sphere-line intersection points
            %and centers of discontinuity and CW
            
            %Point1 = PTSA(1,:)
            d_point1A=dist2points(PTSA(1,:), xyz(i,:));
            d_point1B=dist2points(PTSA(1,:), CWxyz);
            radiusB=CWr;
            
            if ((-d_point1A+radius(i))>-0.0005) && ((-d_point1B+CWr)>-0.0005)
                PTSAB(1,1)=PTSA(1,1);
                PTSAB(1,2)=PTSA(1,2);
                PTSAB(1,3)=PTSA(1,3);
            else
                PTSAB(1,1)=NaN;
                PTSAB(1,2)=NaN;
                PTSAB(1,3)=NaN;
            end
            %Same procedure is done for all other 3 spheres-line intersection points calculated
            %Point 2
            d_point2A=dist2points(PTSA(2,:), xyz(i,:));
            d_point2B=dist2points(PTSA(2,:), CWxyz);
            
            if ((-d_point2A+radius(i))>-0.0005) && ((-d_point2B+radiusB)>-0.0005)
                PTSAB(2,1)=PTSA(2,1);
                PTSAB(2,2)=PTSA(2,2);
                PTSAB(2,3)=PTSA(2,3);
            else
                PTSAB(2,1)=NaN;
                PTSAB(2,2)=NaN;
                PTSAB(2,3)=NaN;
            end
            
            %Point 3
            d_point3A=dist2points(PTSB(1,:), xyz(i,:));
            d_point3B=dist2points(PTSB(1,:), CWxyz);
            if ((-d_point3A+radius(i))>-0.0005) && ((-d_point3B+radiusB)>-0.005)
                PTSAB(3,1)=PTSB(1,1);
                PTSAB(3,2)=PTSB(1,2);
                PTSAB(3,3)=PTSB(1,3);
            else
                PTSAB(3,1)=NaN;
                PTSAB(3,2)=NaN;
                PTSAB(3,3)=NaN;
            end
            %Point 4
            d_point4A=dist2points(PTSB(2,:), xyz(i,:));
            d_point4B=dist2points(PTSB(2,:), CWxyz);
            if ((-d_point4A+radius(i))>-0.0005) && ((-d_point4B+radiusB)>-0.0005)
                PTSAB(4,1)=PTSB(2,1);
                PTSAB(4,2)=PTSB(2,2);
                PTSAB(4,3)=PTSB(2,3);
            else
                PTSAB(4,1)=NaN;
                PTSAB(4,2)=NaN;
                PTSAB(4,3)=NaN;
            end
            
            %--Distances between the points
            d_pp13AB = dist2points(PTSAB(1,:), PTSAB(3,:));
            d_pp14AB = dist2points(PTSAB(1,:), PTSAB(4,:));
            d_pp23AB = dist2points(PTSAB(2,:), PTSAB(3,:));
            d_pp24AB = dist2points(PTSAB(2,:), PTSAB(4,:));
            
            % Erase similar point (save only 1 point for two or more
            % similar points)
            if (d_pp13AB<0.005 && d_pp24AB<0.0005)
                PTSAB(3,1)=NaN;
                PTSAB(3,2)=NaN;
                PTSAB(3,3)=NaN;
                PTSAB(4,1)=NaN;
                PTSAB(4,2)=NaN;
                PTSAB(4,3)=NaN;
                
            end
            if (d_pp23AB<0.0005 && d_pp14AB<0.0005)
                PTSAB(3,1)=NaN;
                PTSAB(3,2)=NaN;
                PTSAB(3,3)=NaN;
                PTSAB(4,1)=NaN;
                PTSAB(4,2)=NaN;
                PTSAB(4,3)=NaN;
                
            end
            
            
            
            
            % Select non-NaN values
            selectAB = ~isnan( PTSAB ) ;
            pointIntersecting=PTSAB(selectAB);
            npointIntersecting= (numel(pointIntersecting))/3;%count number of intersection points
            
            if npointIntersecting>1 %if point int. bigger than 1 a line of intersection is totally included in the
                countCWDisc_inter(i,Set(i))=2;
                
                Point_intersection=[pointIntersecting(1), pointIntersecting(3), pointIntersecting(5), pointIntersecting(2), pointIntersecting(4), pointIntersecting(6)];
                Point_intersCW_disc(i,1:6)=Point_intersection;
                
                X=[pointIntersecting(1); pointIntersecting(2)];Y=[pointIntersecting(3); pointIntersecting(4)];Z=[pointIntersecting(5); pointIntersecting(6)];
                %Mauldon method
                if ((-d_point1B+CWr)>0.0000)%if point are inside the CW is counted as m
                    m_mauldon(i)=1;
                else
                    n_mauldon(i)=1; %else point are counted as n
                end
                if ((-d_point2B+CWr)>-0.0000)
                    m_mauldon(i)=m_mauldon(i)+1;
                else
                    n_mauldon(i)=n_mauldon(i)+1;
                end
                %trace length included in CW
                l(i) = dist2points(Point_intersection(1:3),Point_intersection(4:6));
                l(i) = sqrt((pointIntersecting(1)-pointIntersecting(2))^2 + (pointIntersecting(3)-pointIntersecting(4))^2 + (pointIntersecting(5)-pointIntersecting(6))^2);
                L(i) = dist2points(PTSA(1,:), PTSA(2,:));%ENTIRE TRACE OF DISCONTINUITY
                %---TRACE LENGTH (L) OF ENTIRE DISCONTINUITY (NOT CUTTED FROM CW)
                for j= 1: max(Set)
                    if Set(i) == j
                        Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};% define color for set from 1 to 7
                        drawEdge3d(Point_intersection, 'color', Color{j}, 'linewidth', 4);%plot lines of intersection
                        hold on;grid on
                        
                        % % % % % %                         filename_mod=strrep(fn,'.csv','');
                        % % % % % %                         %DXF names are saved in the same way
                        % % % % % %                         dxf_name=(['Set', num2str(Set(i)),'_Trace',num2str(i),'_',num2str(filename_mod),'.dxf']);
                        % % % % % %                         FID=dxf_open([pn,'results',dxf_name]);
                        % % % % % %                         FID = dxf_set(FID,'Color',Color{j});
                        % % % % % %                         dxf_polyline(FID, X, Y, Z);
                        % % % % % %                         dxf_close(FID);
                    end
                end
            end
            
        else
            
        end
        
    end
    
end
%savefig(fullfile(pathIntersection,'CW_tracemap'))
end