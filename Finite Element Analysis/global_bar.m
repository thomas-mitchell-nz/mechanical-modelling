function [ Khat_e, Lambda_e ] = global_bar ( K_e, alpha )
% GLOBAL_BAR Takes results from LOCAL_BAR.m to generate the transformation
% matrix and perform a mutliplication to get the element stiffness matrix
% in global coordinates 

Lambda_e = [cosd(alpha) sind(alpha) 0 0; 0 0 cosd(alpha) sind(alpha)];
Khat_e = Lambda_e'*K_e*Lambda_e;

end