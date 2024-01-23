function [ Khat_e, Lambda_e ] = global_frame ( K_e, alpha )
% GLOBAL_FRAME Returns element stiffness matrix in global coordinates
% and the corresponding transformation matrix

c = cosd(alpha);
s = sind(alpha);

Lambda_e = [
    c s 0 0 0 0;
    -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
    0 0 0 -s c 0;
    0 0 0 0 0 1
    ];

Khat_e = Lambda_e'*K_e*Lambda_e;

end