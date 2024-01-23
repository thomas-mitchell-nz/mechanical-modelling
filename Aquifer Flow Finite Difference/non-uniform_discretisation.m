clc;close all;clear;

% Parameters
L = 80;
W = 40;
R = 4;

cen_error_last = 0;
max = 99;
min = 9;

data = readmatrix('COMSOL_DATA.txt','headerlines',5);
xFEA = data(:,1);
yFEA = data(:,3);

for f=min:2:max
    Ny = f;
    Nx = Ny;
    str = 0.9;
    nxm=Nx-1;
    tstr3=sinh(str);
    xc(1)=0.0;

    for kc=2:Nx
        z2dp=(2*kc-Nx-1)/(nxm);
        xc(kc)=(1+sinh(str*z2dp)/tstr3)*0.5*L;
    end
        
    nym=Ny-1;
    yc(1)=0.0;
    for kc=2:Ny
        z2dp=(2*kc-Ny-1)/(nym);
        yc(kc)=(1+sinh(str*z2dp)/tstr3)*0.5*W;
    end
      
    [X,Y]=meshgrid(xc,yc);
    
    % Reset parameters
    H = zeros(Ny, Nx);
    rel_error = zeros(Ny, Nx);
    for j=1:Ny
        H(j, 1) = 18;
        H(j, end) = 12;
    end
    K = zeros(Ny, Nx);
    Dx = L/(Nx-1);
    Dy = W/(Ny-1);
    for j=1:(Ny)
        for i=1:(Nx)
            r = sqrt((i-(Nx/2))^2 + (j-(Ny/2))^2);
            if (r < R/(Dy))  
                K(j,i) = 2e-4;
            elseif (r < R/Dx)
                K(j,i) = 2e-4;
            else
                K(j,i) = 2e-3;
            end
        end
    end
    
    
    % Liebmann method
    Hold = H;
    loop = 1;
    iter = 1;
    epsi = 1e-3; % stop6ping criterion for Liebmann Method
    lambda = 1.7; % over-relaxation parameter
    while loop
        for j=1:Ny
            for i=2:(Nx-1)
                
                if (j ==1)
                    K_fx = (K(j,i+1) + K(j,i))/2;
                    K_bx = (K(j,i) + K(j, i-1))/2;
                    K_fy = (K(j+1, i) + K(j,i))/2;
                    K_by = (K(j,i) + K(j+1,i))/2;
    
                    X_cd = (X(j,i+1) - X(j,i));
                    Y_cd = (Y(j+1,i) - Y(j,i));
                
                    X_f = (X(j,i+1) - X(j,i));
                    X_b = (X(j,i) - X(j,i-1));
    
                    Y_f = (Y(j+1,i) - Y(j,i));
                    Y_b = (Y(j,i) - Y(j+1,i));
                    
                elseif (j == Ny)
                    K_fx = (K(j,i+1) + K(j,i))/2;
                    K_bx = (K(j,i) + K(j, i-1))/2;
                    K_fy = (K(j-1, i) + K(j,i))/2;
                    K_by = (K(j,i) + K(j-1,i))/2;
    
                    X_cd = (X(j,i+1) - X(j,i));
                    Y_cd = (Y(j-1,i) - Y(j,i));
                
                    X_f = (X(j,i+1) - X(j,i));
                    X_b = (X(j,i) - X(j,i-1));
    
                    Y_f = (Y(j-1,i) - Y(j,i));
                    Y_b = (Y(j,i) - Y(j-1,i));
                    
                else
                    K_fx = (K(j,i+1) + K(j,i))/2;                
                    K_bx = (K(j,i) + K(j, i-1))/2;
                    K_fy = (K(j+1, i) + K(j,i))/2;
                    K_by = (K(j,i) + K(j-1,i))/2;
    
                    X_cd = (X(j,i+1) - X(j,i));
                    Y_cd = (Y(j+1,i) - Y(j,i));
                
                    X_f = (X(j,i+1) - X(j,i));
                    X_b = (X(j,i) - X(j,i-1));
    
                    Y_f = (Y(j+1,i) - Y(j,i));
                    Y_b = (Y(j,i) - Y(j-1,i));
                    
                end
   
                if (j == 1) % BOTTOM
                    %H(j,i) =  lambda* (((2/X_cd) * ((K_fx/X_f) * H(j, i+1)+(K_bx/X_b) * H(j, i-1)) + (2/Y_cd) * ((K_fy/Y_f)* H(j+1,i) + (K_by/Y_b) * H(j+1, i))) / ((2 * X_cd) * ((K_fx/X_f) * (K_bx/X_b)) + (2/Y_cd)* ((K_fy/Y_f) + (K_by/Y_b)))) +(1-lambda)*Hold(j,i);
                    A = (2 / X_cd) * ((K_fx/X_f) * H(j, i+1) + (K_bx/X_b) * H(j, i-1));
                    
                    B = (2 / Y_cd) * ((K_fy/Y_f) * H(j+1, i) + (K_by/Y_b) * H(j+1, i));
    
                    C = (2 / X_cd) * ((K_fx/X_f) + (K_bx/X_b));
                    D = (2 / Y_cd) * ((K_fy/Y_f) + (K_by/Y_b));
                    
                    H(j,i) =  lambda * (A + B) / (C + D) + (1-lambda) * Hold(j,i);
                elseif (j == Ny) % TOP
                    %H(j,i) =  lambda* (((2/X_cd) * ((K_fx/X_f) * H(j, i+1)+(K_bx/X_b) * H(j, i-1)) + (2/Y_cd) * ((K_fy/Y_f)* H(j-1,i) + (K_by/Y_b) * H(j-1, i))) / ((2 * X_cd) * ((K_fx/X_f) * (K_bx/X_b)) + (2/Y_cd)* ((K_fy/Y_f) + (K_by/Y_b)))) +(1-lambda)*Hold(j,i);
                    A = (2 / X_cd) * ((K_fx/X_f) * H(j, i+1) + (K_bx/X_b) * H(j, i-1));
    
                    B = (2 / Y_cd) * ((K_fy/Y_f) * H(j-1, i) + (K_by/Y_b) * H(j-1, i));
    
    
                    C = (2 / X_cd) * ((K_fx/X_f) + (K_bx/X_b));
                    D = (2 / Y_cd) * ((K_fy/Y_f) + (K_by/Y_b));
                    
                    H(j,i) =  lambda * (A + B) / (C + D) + (1-lambda) * Hold(j,i);
                else
                    %H(j,i) =  lambda* ((  (2/X_cd) * ((K_fx/X_f) * H(j, i+1) + (K_bx/X_b)* H(j, i-1)) + (2/Y_cd) * ((K_fy/Y_f)* H(j+1,i) + (K_by/Y_b) * H(j-1, i)))    / ((2 * X_cd) * ((K_fx/X_f) * (K_bx/X_b)) + (2/Y_cd)* ((K_fy/Y_f) + (K_by/Y_b)))) +(1-lambda)*Hold(j,i);
                    
                   
    
                    A = (2 / X_cd) * ((K_fx/X_f) * H(j, i+1) + (K_bx/X_b) * H(j, i-1));
                    B = (2 / Y_cd) * ((K_fy/Y_f) * H(j+1, i) + (K_by/Y_b) * H(j-1, i));
                    C = (2 / X_cd) * ((K_fx/X_f) + (K_bx/X_b));
                    D = (2 / Y_cd) * ((K_fy/Y_f) + (K_by/Y_b));
                    
                    H(j,i) =  lambda * (A + B) / (C + D) + (1-lambda) * Hold(j,i);
    
                end
                rel_error(j,i) = 100*(H(j,i)-Hold(j,i))/H(j,i);
            end
        end
        norm_error(iter) = norm(rel_error,2);
        if (norm_error(iter)<epsi)||(iter>5000)
            loop = 0;
        else
            iter = iter+1;
            Hold = H;
        end
    end

    
    c_val = H((Ny/2 + 0.5),(Nx/2 + 0.5));
    conv(f/2 - 3.5) = c_val;
    cen_error_last = c_val;
end

conv = conv - c_val;

surf(X,Y,H)
xlabel('x (m)')
ylabel('y (m)')
zlabel('H (m)')

figure
hold on
[C,h] = contour(X,Y,H,10);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
xlabel('x (m)')
ylabel('y (m)')
[qx,qy] = gradient(H,Ny,Nx);
qx = -qx;
qy = -qy;
quiver(X,Y,qx,qy)

figure
x = linspace(min + 1,max,length(conv)-1);
plot(x,abs(conv(2:length(conv))))
xlabel('Number of Nodes Along x')
ylabel('Relative Error (m)')

figure
plot(log10(norm_error),'*r')
xlabel('Iteration number')
ylabel('Log10(Norm Percent Relative Error)')

figure
hold on
plot(xFEA,yFEA, '*')
plot(linspace(1,80,max),H(max/2+0.5,:))
xlabel('Length (m)')
ylabel('Hydraulic Head (m)')
legend(["COMSOL Study","MATLAB study"])
