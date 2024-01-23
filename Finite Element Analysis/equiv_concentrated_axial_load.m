function [ f_CAL,F_CAL,Q_CAL ] = equiv_concentrated_axial_load ( pbar,L,a,Lambda,AM )
% EQUIV_CONCENTRATED_AXIAL_LOAD Returns equivalent nodal loading induced by 
% distributed axial load

f_CAL = pbar*[1-a/L; 0; 0; a/L; 0; 0]; % Equivalent load in local coords

F_CAL = Lambda'*f_CAL; % Global coords 

Q_CAL = AM*F_CAL; % Structural contribution

end