figure;
imagesc(x, z, w);
set(gca, 'YDir', 'normal');
xlabel('$x$ (m)','interpreter','latex','FontSize',16); ylabel('$z$ (m)','interpreter','latex','FontSize',16);
title('Heat map of $w(x,z)$','interpreter','latex','FontSize',16);
colorbar;

figure;
imagesc(x, z, u);
set(gca, 'YDir', 'normal');
xlabel('$x$ (m)','interpreter','latex','FontSize',16); ylabel('$z$ (m)','interpreter','latex','FontSize',16);
title('Heat map of $U_0+u(x,z)$','interpreter','latex','FontSize',16);
colorbar;


figure;
contourf(x, z, w, 20);
set(gca, 'YDir', 'normal');
xlabel('$x$ (m)','interpreter','latex','FontSize',16); ylabel('$z$ (m)','interpreter','latex','FontSize',16);
xlim([2e4,3e4])
ylim([80,460])
title('Contour plot of $w(x,z)$','interpreter','latex','FontSize',16);
colorbar;
