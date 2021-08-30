function [N] = attitude2normal(dipdir,dip,Set)
% This function calculates the normal vector (N) from the attitude
% (dipdirection, dip).

for r=1:length(Set)
    % Unit vector V, which starts at the Origin, in the orientation of the
    % dip/dip direction
    V(:,:,r)=[cosd(90-dipdir(r))*cosd(dip(r)) sind(90-dipdir(r))*cosd(dip(r)) -sind(dip(r))];
    % Unit vector O, which origin is at the origo, which is perpendicular to
    % the dip diprecion and is horizontal
    O(:,:,r)=[-sind(90-dipdir(r)) cosd(90-dipdir(r)) 0];
end
% Normal vector N of the surface defined with the vectors V and O
N=zeros(length(Set),3);

for k=1:length(Set)
    N(k,:)=cross(V(:,:,k),O(:,:,k));
end
end