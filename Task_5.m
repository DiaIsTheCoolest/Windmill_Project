

function dwdz = richardson_dz(w,delta_z,N_x,N_z)
%using Richardson extrapolation to compute dw/dz
    dwdz = zeros(N_z+1,N_x);
    k = 3:N_z-1;
    %centred differences
    d1 = (w(k+1,:)-w(k-1,:))./(2*delta_z);
    d2 = (w(k+2,:)-w(k-2,:))./(4*delta_z);
    %Richardson extrapolation
    dwdz(k,:) = (4*d1-d2)/3;
    %endpoints
    dwdz(2,:) = (w(3,:)-w(1,:))/(2*delta_z);
    dwdz(end-1,:) = (w(end,:)-w(end-2,:))/(2*delta_z);%centred differences
    dwdz(1,:) = (w(2,:)-w(1,:))/delta_z; %forward difference
    dwdz(end,:) = zeros(1,N_x); %initial condition
end

dwdz = richardson_dz(w,delta_z,N_x,N_z);

% Integrate along x to get u
x = linspace(0, L_x, N_x);
u = U0+cumtrapz(x, -dwdz, 2); % integrate -dw/dz along x

u1 = u - U0;

%--------------------------------------------------------------------------

% making contour plots


% w contour plot
figure;
contourf(x, z, w, 20);
set(gca, 'YDir', 'normal');
xlabel('$x$ (m)','interpreter','latex'); 
ylabel('$z$ (m)','interpreter','latex');
title('Contour of $w(x,z)$','interpreter','latex');
colorbar;

% u contour plot
figure;
contourf(x, z, u, 20);
set(gca, 'YDir', 'normal');
xlabel('$x$ (m)','interpreter','latex'); 
ylabel('$z$ (m)','interpreter','latex');
title('Contour of $U_0+u(x,z)$','interpreter','latex');
colorbar;

% vector field
% (not used in the report in the end - the info seen in this can be seen through the heat maps)
figure;
skip = 2;
quiver(X(1:skip:end,1:skip:end), Z(1:skip:end,1:skip:end), u1(1:skip:end,1:skip:end), w(1:skip:end,1:skip:end), LineWidth=0.5, ShowArrowHead="on", MaxHeadSize=0.0005);
% quiver(X, Z, u1, w, LineWidth=0.5, ShowArrowHead="on", MaxHeadSize=0.0001);
xlabel('$x$ (m)', 'interpreter', 'latex');
ylabel('$z$ (m)', 'interpreter', 'latex');
title('Velocity vector field of $u$ and $w$, without $U0$', 'interpreter', 'latex');
axis tight
set(gca, 'FontSize', 12);
