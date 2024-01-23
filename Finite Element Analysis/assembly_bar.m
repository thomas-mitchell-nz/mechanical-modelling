function [ K_G ] = assembly_bar ( Khat_e, A_e )
% ASSEMBLY_BAR Returns the element's contribution to the global stiffness
% matrix

K_G = A_e*Khat_e*A_e';

end