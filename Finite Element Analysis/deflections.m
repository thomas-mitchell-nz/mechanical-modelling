function [ u,v,D_XG,D_YG ] = deflections ( d_e,L,x_e,alpha )
% DEFLECTIONS Returns the axial and transverse deflections of an element

% Transverse shape functions
N1 = 1-3*x_e.^2/L^2+2*x_e.^3/L^3;
N2 = x_e.^3/L^2-2*x_e.^2/L+x_e;
N3 = 3*x_e.^2/L^2-2*x_e.^3/L^3;
N4 = x_e.^3/L^2-x_e.^2/L;

% Element transverse displacement
v = N1*d_e(2)+N2*d_e(3)+N3*d_e(5)+N4*d_e(6);

% Axial shape functions
psi1 = (1-x_e/L);
psi2 = x_e/L;

% Element axial displacement 
u = psi1*d_e(1) + psi2*d_e(4);

% Transform to global coords
D_YG = u*sind(alpha)+v*cosd(alpha);
D_XG = u*cosd(alpha)-v*sind(alpha);

end