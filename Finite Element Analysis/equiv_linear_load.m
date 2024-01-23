function [ f_LVL,F_LVL,Q_LVL ] = equiv_linear_load ( wbar,L,Lambda,AM )
% EQUIV_LINEAR_LOAD Returns equivalent nodal loading induced by LVL

% Equivalent load in local coords
f_LVL = [0;
         3*wbar*L/20;
         wbar*L^2/30;
         0;
         7*wbar*L/20;
         -wbar*L^2/20];

F_LVL = Lambda'*f_LVL; % Global coords

Q_LVL = AM*F_LVL; % Structural contribution

end