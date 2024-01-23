function [ K_e ] = local_bar ( E,A,L )
% LOCAL_BAR Given E,A,L generate the element stiffness matrix in local
% cooridnates 

K_e = E*A/L*[1 -1; -1 1];

end