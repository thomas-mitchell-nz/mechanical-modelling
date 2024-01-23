function [ Deflected_XG, Deflected_YG ] = Plot_deflected_bar ...
    (node1XG,node1YG,node2XG,node2YG,D_e,Disp_mag)
% PLOT_DEFLECTED_BAR Generate a plot for the deflected shape of a generic
% bar element

% Transfrom displacements from local to global coords
Deflections_XG = [D_e(1) D_e(3)];
Deflections_YG = [D_e(2) D_e(4)];

% Undeflected baseline of the element
Undeflected_baseline_XG = linspace(node1XG,node2XG,2);
Undeflected_baseline_YG = linspace(node1YG,node2YG,2);

% Deflected position of all points within the element
Deflected_XG = Undeflected_baseline_XG+Disp_mag*Deflections_XG;
Deflected_YG = Undeflected_baseline_YG+Disp_mag*Deflections_YG;

% Plot
plot(Undeflected_baseline_XG,Undeflected_baseline_YG,'b.-')
hold on
plot(Deflected_XG,Deflected_YG,'r.-')

end