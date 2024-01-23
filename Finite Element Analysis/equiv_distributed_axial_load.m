function [ f_DAL,F_DAL,Q_DAL ] = equiv_distributed_axial_load ( pbar,L,Lambda,AM )
% EQUIV_DISTRIBUTED_AXIAL_LOAD Returns equivalent nodal loading induced by 
% distributed axial load

f_DAL = pbar*[L/2; 0; 0; L/2; 0; 0]; % Equivalent load in local coords

F_DAL = Lambda'*f_DAL; % Global coords

Q_DAL = AM*F_DAL; % Structural contribution

end