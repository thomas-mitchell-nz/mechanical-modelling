function [ f_PL,F_PL,Q_PL ] = equiv_point_load ( wbar,L,a,Lambda,AM )
% EQUIV_POINT_LOAD Returns equivalent nodal loading induced by point load

% Equivalent load in local coords
f_PL = [0;
        1-3*a^2/L^2+2*a^3/L^3;
        a^3/L^2-2*a^2/L+a;
        0;
        3*a^2/L^2-2*a^3/L^3;
        a^3/L^2-a^2/L]*wbar;

F_PL = Lambda'*f_PL; % Global coords

Q_PL = AM*F_PL; % Structural contribution

end