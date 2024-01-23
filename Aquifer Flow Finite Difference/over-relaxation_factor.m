clc; 
close all;
clear;

% Parameters
L = 80; % m
W = 40; % m
dx = 5; % grid point spacing in x 
dy = 5; % grid point spacing in y
Nx = L/dx+1; % num of grid points in x
Ny = W/dy+1; % num of grid points in y
k = 0.02; % thermal conductivity

% Mesh grid
x = linspace(0,L,Nx); % x-coordinate of grid points (m)
y = linspace(0,W,Ny); % y-coordinate of grid points (m)
[X,Y] = meshgrid(x,y); % matrix of x and y coordinates (m), (useful for plotting)

% Preallocate arrays
H = zeros(Ny, Nx);
rel_error = zeros(Ny, Nx);
lambdas = zeros(10,1);
num_iterations = zeros(10,1);

% Boundary conditions
for j=1:Ny
    H(j, 1) = 18;
    H(j, end) = 12;
end

% Liebmann method
Hold = H;
loop = 1;
iter = 1;
epsi = 1e-3; % stopping criterion for Liebmann Method

k = 1; % iterator
for lambda=1:0.01:2 % over-relaxation parameters
    while loop
        for j=1:(Ny)
            for i=2:(Nx-1)
                if (j == 1)
                    H(j,i) =  lambda*(0.25*(H(j,i+1)+H(j,i-1)+ 2*H(j+1,i)))+(1-lambda)*Hold(j,i);
                elseif (j == Ny)
                    H(j,i) =  lambda*(0.25*(H(j,i+1)+H(j,i-1)+ 2*H(j-1,i)))+(1-lambda)*Hold(j,i);
                else
                    H(j,i) =  lambda*(0.25*(H(j,i+1)+H(j,i-1)+H(j+1,i)+H(j-1,i)))+(1-lambda)*Hold(j,i);
                end
                rel_error(j,i) = 100*(H(j,i)-Hold(j,i))/H(j,i);
            end
        end
        norm_error(iter) = norm(rel_error,2);
        if (norm_error(iter)<epsi)||(iter>2000)
            loop = 0;
        else
            iter = iter+1;
            Hold = H;
        end
    end

    % Table data
    num_iterations(k) = iter;
    lambdas(k) = lambda;
    k = k+1;
    
    % Reset parameters
    Hold = H;
    loop = 1;
    iter = 1;
    H = zeros(Ny, Nx); % pressure
    rel_error = zeros(Ny, Nx);
    for j=1:Ny
        H(j, 1) = 18;
        H(j, end) = 12;
    end
end

% Table
relax_table = table;
relax_table.Relaxation_factor = lambdas;
relax_table.Number_of_iterations = num_iterations;
disp(relax_table)

% Optimal lambda
[M,I]=min(num_iterations);
disp('Optimal lambda value = ');disp(lambdas(I));