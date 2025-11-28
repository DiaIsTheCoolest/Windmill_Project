t_start = 0;
t_end = 3600;
delta_t = 1;
N = floor((t_end - t_start) / delta_t);
timeVector = t_start:delta_t:t_end;

X_trajectory = zeros(1, N+1);
Z_trajectory = zeros(1, N+1);



U = u + U0;
W = w;


function velocity = get_velocities(MATRIX,x,y, X,Z)
     velocity = interp2(X, Z, MATRIX, x, y, "linear", 0); %Interpolates an estimate for the velocity, if the particle does not sit exactly on a grid point
end



X0 = 280;
Z0 =100;

X_current = X0;
Z_current = Z0;

for i = 1:N
    X_trajectory(i) = X_current;
    Z_trajectory(i) = Z_current;




    u_n = get_velocities(U, X_current, Z_current, X,Z);
    w_n = get_velocities(W, X_current, Z_current, X, Z);

    X_star_new = X_current + u_n*delta_t;
    Z_star_new = Z_current + w_n*delta_t;
    
    
    u_star_new = get_velocities(U,X_star_new, Z_star_new, X, Z);
    w_star_new = get_velocities(W,X_star_new, Z_star_new, X, Z);


    X_new = X_current + delta_t * (u_n + u_star_new)/2;
    Z_new = Z_current + delta_t * (w_n+w_star_new)/2;
        

    X_current = X_new;
    Z_current = Z_new;
end
X_trajectory(N) = X_current;
Z_trajectory(N) = Z_current;


figure
plot(timeVector, Z_trajectory)
xlabel('Time (s)')
ylabel('Z (m)')
title('Z Trajectory vs Time')
grid on


xlim([1,45])

figure
plot(timeVector, X_trajectory)
xlabel('Time (s)')
ylabel('X (m)')
title('X Trajectory vs Time')
grid on


xlim([1,45])



figure
plot(X_trajectory, Z_trajectory)
hold on
quiver(x,z,U, W)
xlabel('X (m)')
ylabel('Z (m)')
title('X Trajectory vs Z trajectory')
grid on





