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
lambda = 1.5; % over-relaxation factor

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
    c_error(iter) = H(5,9);
    if (norm_error(iter)<epsi)||(iter>1000)
        loop = 0;
    else
        iter = iter+1;
        Hold = H;
    end
end
disp('Number of iteration to convergence = ');disp(iter)

surf(X,Y,H)
xlabel('x (m)')
ylabel('y (m)')
zlabel('H (m)')

figure
[C,h] = contour(X,Y,H,10);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
xlabel('x (m)')
ylabel('y (m)')

% Solution gradient on contour
[qx,qy] = gradient(H,dx,dx);
qx = -k*qx;
qy = -k*qy;
hold on
quiver(x,y,qx,qy)

figure
plot(log10(norm_error),'*r')
xlabel('Iteration number')
ylabel('Log10(Norm Percent Relative Error)')

figure
plot(c_error)
yline(15)
xlabel('Number of Iterations')
ylabel('Hydraulic Head (m)')
