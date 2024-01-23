function [ K_G ] = assembly_frame ( Khat_e, A_e )
% ASSEMBLY_FRAME Returns the element's contribution to the global stiffness
% matrix

K_G = A_e*Khat_e*A_e';

end