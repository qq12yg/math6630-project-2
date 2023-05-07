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

% Crank-Nicolson method
U = zeros(N, M+1);
U(:,1) = u0;

A = eye(N);
B = eye(N);

for j = 1:N
    for k = 1:N
        A(j,k) = A(j,k) - dt/4 * (sin(2*pi*x(j)) * (2*pi*i*(k-1)) + 0.5 * (-4*pi^2*(k-1)^2));
        B(j,k) = B(j,k) + dt/4 * (sin(2*pi*x(j)) * (2*pi*i*(k-1)) + 0.5 * (-4*pi^2*(k-1)^2));
    end
end

for m = 1:M
    U(:,m+1) = A \ (B * U(:,m));
end

% Plot the solution
surf(t, x, real(U));
xlabel('Time');
ylabel('Space');
zlabel('u(x, t)');
title('Crank-Nicolson solution');