function [ K_e ] = local_frame ( E,I,A,L )
% LOCAL_FRAME Given E,I,A,L generate the element stiffness matrix in local
% cooridnates 

B = A*L^2/I;

K_e = E*I/L^3 * [
    B 0 0 -B 0 0;
    0 12 6*L 0 -12 6*L;
    0 6*L 4*L^2 0 -6*L 2*L^2;
    -B 0 0 B 0 0;
    0 -12 -6*L 0 12 -6*L;
    0 6*L 2*L^2 0 -6*L 4*L^2
    ];

end