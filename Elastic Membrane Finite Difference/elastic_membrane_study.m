% ENME302
% Thomas Mitchell
% 1D FEA study of elastic membrane with linear and quadratic element models

clear
close all
clc

% parameters
W = 0.01; % m
f = 100; % N/m3
tau0 = 1; % Pa

% data from COMSOL 1D study
quad = readmatrix('1d_quad_plot.txt','headerlines',9);
lin = readmatrix('1d_linear_plot.txt','headerlines',9);
xc_quad = quad(:,1);
ugc_quad = quad(:,2);
xc_lin = lin(:,1);
ugc_lin = lin(:,2);


%% Linear elements

% Linear shape element
syms x x1 x2 u1 u2 a0 a1 % define the symbolic variables we are using
u = a0 + a1*x; % interpolation function (linear)
tau = 1 + 0.9*sin(sqrt(2)*pi*x/W); % stiffness equation
eqn1 = subs(u,x,x1) == u1; % evaluating at left node
eqn2 = subs(u,x,x2) == u2; % evaluating at right node
cons = solve([eqn1 eqn2], [a0 a1]); % solve simultaneous equations for constants a and b
u = subs(u,[a0 a1], [cons.a0 cons.a1]); % substitute/plug into interpolation function
[N,~] = coeffs(u, [u1 u2]); % collect coefficients for defining shape functions
N1 = N(1); % N1 = (x2-x)/(x2-x1)
N2 = N(2); % N2 = (x-x1)/(x2-x1)

% Element equations
Ke = sym(zeros(2)); % element stiffness matrix
Fe = sym(zeros(2,1)); % element rhs forcing vector
for i = 1:2 % 2 dof
    eleq = int(tau*diff(u,x)*diff(N(i),x),x,x1,x2)-int(f*N(i),x,x1,x2); % element equation
    [coef,~] = coeffs(eleq, [u1 u2]); % collect coefficients
    Ke(i,:) = coef(1:2); % int(diff(u,x)*diff(N(i),x),x,x1,x2)
    Fe(i)= -coef(3); % int(f*N(i),x,x1,x2)
end

% Global equations 
noe = 2; % number of elements
dof = noe+1; % number of degrees of freedom
x = linspace(0,W,dof); % x-coordinates of dof
Kg = zeros(dof); % global stiffness matrix
Fg = zeros(dof,1); % global rhs forcing vector
for i=1:noe % assemble each set of element equations
    Kg(i:i+1,i:i+1) = Kg(i:i+1,i:i+1) + subs(Ke,[x1 x2], [x(i) x(i+1)]); % global contribution
    Fg(i:i+1) = Fg(i:i+1) + subs(Fe,[x1 x2], [x(i) x(i+1)]) ; % global contribution
end 

% Boundary conditions 
ua = 0; % Dirichlet boundary condition along a 
Fg(1) = Fg(1) - Kg(1,1)*ua;
Kg(1,1) = 1;
Fg(2) = Fg(2) - Kg(2,1)*ua;
Kg(2,1) = 0;
ub = 0; % Dirichlet boundary condition along b
Fg(dof) = Fg(dof) - Kg(dof,dof)*ub;
Kg(dof,dof) = -1*(1+0.9*sin((sqrt(2)*pi*W)/W));
Fg(dof-1) = Fg(dof-1) - Kg(dof-1,dof)*ub;
Kg(dof-1,dof) = 0;

% Solution
x_lin = x; % independent variable
kg_lin = Kg; % global stiffness matrix
ubc = Kg\Fg; % distribution with BCs applied
ug_lin = [ua ubc(2:dof-1)' ub]; % global vector

%% Quadratic elements

% Quadratic shape element
syms x x1 x2 x3 u1 u2 u3 a b c % define the symbolic variables we are using
u = a*x^2 + b*x + c; % interpolation function (quadratic)
tau = 1 + 0.9*sin(sqrt(2)*pi*x/W); % stiffness equation
eqn1 = subs(u,x,x1) == u1; % evaluating at node 1
eqn2 = subs(u,x,x2) == u2; % evaluating at node 2
eqn3 = subs(u,x,x3) == u3; % evaluating at node 3
cons = solve([eqn1 eqn2 eqn3], [a b c]); % solve simultaneous equations for constants a, b and c
u = subs(u,[a b c], [cons.a cons.b cons.c]); % substitute/plug into interpolation function
[N,~] = coeffs(u, [u1 u2 u3]); % collect coefficients for defining shape functions
N1=N(1); % N1 = 2*xi^2 - 3*xi + 1
N2=N(2); % N2 = -4*xi^2 + 4*xi
N3=N(3); % N3 = 2*xi^2 - xi

% Element equations
Ke = sym(zeros(3)); % element stiffness matrix
Fe = sym(zeros(3,1)); % element rhs forcing vector
for i = 1:3 % 3 dof
    eleq = int(tau*diff(u,x)*diff(N(i),x),x,x1,x3)-int(f*N(i),x,x1,x3); % element equation
    [coef,~] = coeffs(eleq, [u1 u2 u3]); % collect coefficients
    Ke(i,:) = coef(1:3); % int(diff(u,x)*diff(N(i),x),x,x1,x2)
    Fe(i)= -coef(4); % int(f*N(i),x,x1,x2)
end

% Global equations 
noe = 2; % number of elements
dof = 2*noe+1; % number of degrees of freedom
x = linspace(0,W,dof); % x-coordinates of dof
Kg = sym(zeros(dof)); % global stiffness matrix
Fg = sym(zeros(dof,1)); % global rhs forcing vector
for i=1:2:2*noe % assemble each set of element equations
   Kg(i:i+2,i:i+2) = Kg(i:i+2,i:i+2) + subs(Ke,[x1 x2 x3], [x(i) x(i+1) x(i+2)]); % global contribution
   Fg(i:i+2) = Fg(i:i+2) + subs(Fe,[x1 x2 x3], [x(i) x(i+1) x(i+2)]);  % global contribution
end 

% Boundary conditions 
ua = 0; % Dirichlet boundary condition along a 
%Fg(1) = Fg(1) - Kg(1,1)*ua;
Kg(1,1) = 1;
Kg(2,1) = 0;
Kg(3,1) = 0;
%Fg(2) = Fg(2) - Kg(2,1)*ua;
Kg(dof-2,dof) = 0;
ub = 0; % Dirichlet boundary condition along b
%Fg(dof) = Fg(dof) - Kg(dof,dof)*ub;
Kg(dof,dof) = -1*(1+0.9*sin((sqrt(2)*pi*W)/W));
%Fg(dof-1) = Fg(dof-1) - Kg(dof-1,dof)*ub;
Kg(dof-1,dof) = 0;

% Solution
x_quad = x; % independent variable
Kg_quad = Kg; % global stiffness matrix
ubc = Kg\Fg; % distribution with BCs applied
ug_quad = [ua ubc(2:dof-1)' ub]; % global vector

%% Plot

% Quadratic fit
el1x = double(x_quad(1:3)');
el1y = double(ug_quad(1:3)');
el2x = double(x_quad(3:5)');
el2y = double(ug_quad(3:5)');
[el1,~] = fit(el1x,el1y,'poly2');
[el2,~] = fit(el2x,el2y,'poly2');

% Figure
hold on
pl_lin = plot(x_lin,ug_lin);
pl_quad = plot(el1,el1x,el1y); % element 1 of quadratic solution
plot(el2,el2x,el2y); % element 2 of quadratic solution
plc_quad = plot(xc_quad,ugc_quad,'*');
plc_lin = plot(xc_lin,ugc_lin,'*');

xlabel('x (m)'); 
ylabel('u (m)');
legend([pl_lin plc_lin pl_quad(2) plc_quad],["Linear MATLAB" "Linear COMSOL" "Quadratic MATLAB" "Quadratic COMSOL"])