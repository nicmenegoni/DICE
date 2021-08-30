function [p21, muL2D, stdL2D, maxL2D,minL2D,muL3D, stdL3D, maxL3D,minL3D] = P21calc(CWxyz, CWNxyz, CWr, xyz, Nxyz, radius)
%% INFORMATION
% This function count verify if a discontinuity intersects the scan plane
% of the RoKA code.
% It is based onto another script of the DIEapp suite that perform the P21
%
%--------------INPUT parameters:------------------------------------------
%   CWxyz = center of the Circular Window;
%   CWNxyz = normal vector of the Circular Window;
%   CWr = radius of the Circular Window;
%   xyz = center coordinates of the disconutinuity;
%   Nxyz = normal of the discontinuity plane;
%   radius = radius of the disconunuity disc;
%
%-------------OUTPUT parameters:------------------------------------------
%   countCWDisc_inter = vector array where the results of the intersection
%   test is stored (0 = no intersection, 1 = intersection)
%
%% Circular window - discontinuities intersections
nplane = length(radius);
% Preallocating memory
lineAB=zeros(nplane,6);
Int_Length=zeros(nplane,1);

for ii=1:nplane
    %%Compute intersection between two plane
    tol = 1e-14;% setting the angle cutoff for which pseudo parallel planes (parallel at the CW plane) are not consider.
    if abs(sum(cross(Nxyz(ii,:), CWNxyz, 2))/3) > tol
        % The 10 rows below  are based onto a code achived from the website
        % tbirdal.blogspot.com and are based onto a function created
        % according the the solution prioposed by Krumm (2010).
        %References:
        % Krumm J. (2010). Interesection of Two Planes. Microsoft Reseacrh. Available
        % at: https://www.microsoft.com/en-us/research/publication/intersection-of-two-planes/
        M = [2 0 0 CWNxyz(1) Nxyz(ii,1)
            0 2 0 CWNxyz(2) Nxyz(ii,2)
            0 0 2 CWNxyz(3) Nxyz(ii,3)
            CWNxyz(1) CWNxyz(2) CWNxyz(3) 0 0
            Nxyz(ii,1) Nxyz(ii,2) Nxyz(ii,3) 0 0];
        
        b4 = CWxyz(1).*CWNxyz(1) + CWxyz(2).*CWNxyz(2) + CWxyz(3).*CWNxyz(3);
        b5 = xyz(ii,1).*Nxyz(ii,1) + xyz(ii,2).*Nxyz(ii,2) + xyz(ii,3).*Nxyz(ii,3);
        b = [2*CWxyz(1) ; 2*CWxyz(2) ; 2*CWxyz(3); b4 ; b5];
        
        x = M\b;
        Pline = x(1:3)';%position of the vector line
        Vline = cross(CWNxyz, Nxyz(ii,:));%direction vector of the intersection line
        
        lineAB(ii, 1:6)=[Pline,Vline];%Infinite line representing intersection between two infinte planes of discontinuity and CW
        
        point1AB(ii,:)= xyz(ii,:);% Center of discontinuity
        point2AB= CWxyz;% center of Circular Window
        
        %Distance between disconintuity center and previous intersection line
        %dpointL_PA(ii) = abs(cross((xyz(ii,:)-Pline),(xyz(ii,:)-(Pline+Vline)))) / abs(Pline - (Pline+Vline));
        dpointL_PA(ii) = norm(cross(Vline,(Pline-xyz(ii,:))) /norm(Vline));
        
        %Distance between circular window center and previous intersection line
        dpointL_PB(ii) = norm(cross(Vline,(Pline-CWxyz))/norm(Vline));
        
        if abs(dpointL_PA(ii))<radius(ii) && abs(dpointL_PB(ii))<CWr %if distance is lower than radius of 2 discs (discont. and CW) intersection is considered
            selectAB=zeros(4,4);
            
            %Define Sphere (fitted for the 2 considered discs) for Sphere-Line Intersection
            SPHEREA = [xyz(ii,1), xyz(ii,2), xyz(ii,3),  radius(ii)];%definisco la sfera per la discontinuità
            SPHEREB = [CWxyz(1,1), CWxyz(1,2), CWxyz(1,3),  CWr];%definisco la sfera per la CW
            
            %Use Spehere-Line intersection function of geom3d package
            %(INRA) to find points of intersection
            PTSA= intersectLineSphere(lineAB(ii,:), SPHEREA);% for line and discontinuity sphere
            P_interA(ii,1:6)=[PTSA(1,1:3), PTSA(2,1:3)];
            PTSB = intersectLineSphere(lineAB(ii,:), SPHEREB);% for line and CW sphere
            P_interB(ii,1:6)=[PTSB(1,1:3), PTSB(2,1:3)];
            
            % Calculation of l, the cutted trace length, trace visble
            % inside the CW
            
            %calculation distance between sphere-line intersection points
            %and centers of discontinuity and CW
            
            %Point1 = PTSA(1,:)
            d_point1A = sqrt( (xyz(ii,1)-PTSA(1,1))^2 + (xyz(ii,2)-PTSA(1,2))^2 + (xyz(ii,3)-PTSA(1,3))^2);
            d_point1B = sqrt( (CWxyz(1)-PTSA(1,1))^2 + (CWxyz(2)-PTSA(1,2))^2 + (CWxyz(3)-PTSA(1,3))^2);
            radiusB=CWr;
            
            if ((-d_point1A+radius(ii))>-0.0005) && ((-d_point1B+CWr)>-0.0005)
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
            d_point2A = sqrt( (xyz(ii,1)-PTSA(2,1))^2 + (xyz(ii,2)-PTSA(2,2))^2 + (xyz(ii,3)-PTSA(2,3))^2);
            d_point2B = sqrt( (CWxyz(1)-PTSA(2,1))^2 + (CWxyz(2)-PTSA(2,2))^2 + (CWxyz(3)-PTSA(2,3))^2);
            
            if ((-d_point2A+radius(ii))>-0.0005) && ((-d_point2B+radiusB)>-0.0005)
                PTSAB(2,1)=PTSA(2,1);
                PTSAB(2,2)=PTSA(2,2);
                PTSAB(2,3)=PTSA(2,3);
            else
                PTSAB(2,1)=NaN;
                PTSAB(2,2)=NaN;
                PTSAB(2,3)=NaN;
            end
            
            %Point 3
            d_point3A = sqrt( (xyz(ii,1)-PTSB(1,1))^2 + (xyz(ii,2)-PTSB(1,2))^2 + (xyz(ii,3)-PTSB(1,3))^2);
            d_point3B = sqrt( (CWxyz(1)-PTSB(1,1))^2 + (CWxyz(2)-PTSB(1,2))^2 + (CWxyz(3)-PTSB(1,3))^2);
            
            if ((-d_point3A+radius(ii))>-0.0005) && ((-d_point3B+radiusB)>-0.005)
                PTSAB(3,1)=PTSB(1,1);
                PTSAB(3,2)=PTSB(1,2);
                PTSAB(3,3)=PTSB(1,3);
            else
                PTSAB(3,1)=NaN;
                PTSAB(3,2)=NaN;
                PTSAB(3,3)=NaN;
            end
            %Point 4
            d_point4A = sqrt( (xyz(ii,1)-PTSB(2,1))^2 + (xyz(ii,2)-PTSB(2,2))^2 + (xyz(ii,3)-PTSB(2,3))^2);
            d_point4B = sqrt( (CWxyz(1)-PTSB(2,1))^2 + (CWxyz(2)-PTSB(2,2))^2 + (CWxyz(3)-PTSB(2,3))^2);
            
            if ((-d_point4A+radius(ii))>-0.0005) && ((-d_point4B+radiusB)>-0.0005)
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
                if pointIntersecting(5)>pointIntersecting(6)
                    %set the point with higher Z coord. as Point1 (not Point ii !!)
                    %in the final table Point ii could be not equal to
                    %Point1. Point1 has to be the point with higher Z
                    %to obtain the correct plunge and trend values
                    X=[pointIntersecting(1); pointIntersecting(2)];
                    Y=[pointIntersecting(3); pointIntersecting(4)];
                    Z=[pointIntersecting(5); pointIntersecting(6)];
                else
                    X=[pointIntersecting(2); pointIntersecting(1)];
                    Y=[pointIntersecting(4); pointIntersecting(3)];
                    Z=[pointIntersecting(6); pointIntersecting(5)];
                end
                
                %Calculation of length of the intersection between
                %discontinuity disc and scan-plane
                Int_Length(ii)= sqrt((X(2,1) - X(1,1))^2 + (Y(2,1) - Y(1,1))^2 + (Z(2,1) - Z(1,1))^2);
                
                
                
            end
            
            
        end
        
    end
    
end
muL2D=mean(Int_Length);
stdL2D=std(Int_Length);
maxL2D=max(Int_Length);
minL2D=min(Int_Length);
if max(Int_Length)>0
    muL3D=mean(radius(Int_Length>0));
    stdL3D=std(radius(Int_Length>0));
    maxL3D=max(radius(Int_Length>0));
    minL3D=min(radius(Int_Length>0));
else
    muL3D=0;
    stdL3D=0;
    maxL3D=0;
    minL3D=0;
end
p21=sum(Int_Length)/(pi * CWr^2);%ratio of sum of the traces lengths (inside the circular window)/area of the scan circle
end