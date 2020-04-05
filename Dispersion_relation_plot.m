%% Dispersion relation 3D plot
[X,Y] = meshgrid(-2*pi:0.2:2*pi);
Z = -2 * (cos(X) + 2 * cos(0.5*X).*cos(sqrt(3)/2*Y));
figure(1);
surfc(X,Y,Z);
title("Dispersion relation \epsilon_k", 'Fontsize', 16);
xlabel("k_x", 'Fontsize', 14);
ylabel("k_y", 'Fontsize', 14);

contour(X,Y,Z);
title("Dispersion relation \epsilon_k: contour plot", 'Fontsize', 16);
xlabel("k_x", 'Fontsize', 14);
ylabel("k_y", 'Fontsize', 14);