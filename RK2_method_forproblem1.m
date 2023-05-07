% Parameters
L = 1;
T = 1;
N = 50; % spatial discretization
M = 1000; % time discretization
dx = L/N;
dt = T/M;

% Spatial discretization
x = (0:dx:L-dx).';

% Time discretization
t = 0:dt:T;

% Initial condition
u0 = sin(2*pi*x);

% RK2 method
U = zeros(N, M+1);
U(:,1) = u0;

% Define function for the PDE
F = @(u) -sin(2*pi*x).*gradient(u,dx) + 0.5*gradient(gradient(u,dx),dx);

for m = 1:M
    k1 = F(U(:,m));
    k2 = F(U(:,m) + dt*k1);
    U(:,m+1) = U(:,m) + dt/2 * (k1 + k2);
end

% Plot the solution
surf(t, x, real(U));
xlabel('Time');
ylabel('Space');
zlabel('u(x, t)');
title('RK2 solution');