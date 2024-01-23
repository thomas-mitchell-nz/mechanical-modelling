function [ f_UDL,F_UDL,Q_UDL ] = equiv_uniform_load ( wbar,L,Lambda,AM )
% EQUIV_UNIFORM_LOAD Returns equivalent nodal loading induced by UDL

% Equivalent load in local coords
f_UDL = [0;
         wbar*L/2;
         wbar*L^2/12;
         0;
         wbar*L/2;
         -wbar*L^2/12];

F_UDL = Lambda'*f_UDL; % Global coords

Q_UDL = AM*F_UDL; % Structural contribution

end