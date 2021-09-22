function [p, n] = intersect_planes(p1, n1, p2, n2)
% This funtion is achived from the website tbirdal.blogspot.com and is a
% funtion created according the the solution prioposed by Krumm (2010).
%References:
% Krumm J. (2010). Interesection of Two Planes. Microsoft Reseacrh. Available
% at: https://www.microsoft.com/en-us/research/publication/intersection-of-two-planes/
M = [2 0 0 n1(1) n2(1)
     0 2 0 n1(2) n2(2)
     0 0 2 n1(3) n2(3)
     n1(1) n1(2) n1(3) 0 0
     n2(1) n2(2) n2(3) 0 0
    ];

b4 = p1(1).*n1(1) + p1(2).*n1(2) + p1(3).*n1(3);
b5 = p2(1).*n2(1) + p2(2).*n2(2) + p2(3).*n2(3);
b = [2*p1(1) ; 2*p1(2) ; 2*p1(3); b4 ; b5];

x = M\b;
p = x(1:3)';
n = cross(n1, n2);
end