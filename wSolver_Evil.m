% --------------------------
% Parameters
L_x = 1000;     % domain length in x
L_z = 800;     % domain length in z
x_0 = L_x/4;   % Gaussian center x
z_0 = 150;     % Gaussian center z
A = 0.5;       % forcing amplitude
U0 = 20;       % background velocity
sig_x = 50;    % Gaussian width x
sig_z = 40;    % Gaussian width z
N_x = 128;     % number of x points (FFT-friendly)
N_z = 100;     % number of z intervals

% --------------------------
% Grids
x = linspace(0,L_x,N_x);      % x-grid
z = linspace(0,L_z,N_z+1);    % z-grid

[X,Z] = meshgrid(x,z);

% --------------------------
% Forcing: Gaussian
f = -A .* exp(-(X-x_0).^2/(2*sig_x^2)) .* exp(-(Z-z_0).^2/(2*sig_z^2));

% --------------------------
% FFT along x (columns)
f_hat = fft(f,[],2);   % FFT along x, size (N_z+1 x N_x)

% --------------------------
% Solve ODE along z for each Fourier mode
w_hat = zeros(size(f_hat));  % preallocate

delta_z = L_z / N_z;        % grid spacing in z
N = N_z;

% Finite difference matrix for z-derivatives
a = U0 * (-2/delta_z^2 - (2*pi/L_x)^2); % will overwrite k below
b = U0 / delta_z^2;
main_diag = zeros(N,1);
off_diag  = b * ones(N,1);

for j = 1:N_x
    % Fourier wavenumber
    k = 2*pi*(j-1)/L_x;  % MATLAB indexing: 0..N_x-1
    a = U0 * (-2/delta_z^2 - k^2);
    main_diag(:) = a;
    A = spdiags([off_diag main_diag off_diag], -1:1, N, N);
    % Neumann BC at top
    A(N,N-1) = 2*b;
    
    % Compute df/dz using centered differences
    f_col = f_hat(:,j);   % size N+1
    F_k = zeros(N,1);
    F_k(2:N-1) = -(f_col(3:N) - f_col(1:N-2)) / (2*delta_z);
    F_k(1) = -(f_col(2)-f_col(1))/delta_z;
    F_k(N) = -(f_col(N+1)-f_col(N))/delta_z;
    
    % Solve system
    w_k = A \ F_k;
    w_hat(:,j) = [0; w_k];  % include w(0)=0
end

% --------------------------
% Convert back to real space
w = ifft(w_hat,[],2,'symmetric');  % w(x,z), size (N_z+1 x N_x)

% --------------------------
% Plot solution
figure;
imagesc(x,z,w);
set(gca,'YDir','normal');
xlabel('x'); ylabel('z');
title('w(z,x) from pseudo-spectral solution');
colorbar;
% Plot the contour of the solution
figure;
contourf(x, z, w, 20);
set(gca, 'YDir', 'normal');
xlabel('x'); ylabel('z');
title('Contour of w(z,x)');
colorbar;

