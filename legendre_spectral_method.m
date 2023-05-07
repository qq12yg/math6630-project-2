% Parameters
N = 64; % Number of Legendre polynomials
T = 1; % Final time
dt = 0.001; % Time step
mu = 0.1; % Viscosity

% Spatial grid
x = cos(pi*(0:N)/N)';
dx = diff(x);
x_mid = (x(1:end-1)+x(2:end))/2;
J = diag(2./dx) - diag(1./dx(1:end-1),1) - diag(1./dx(1:end-1),-1);
J = J(2:end-1,2:end-1);
L = pi*diag([0:N])/(2*N);

% Initial condition
u0 = sin(pi*x);

% Legendre-Galerkin solver
A = eye(N) - 0.5*dt*mu*L^2;
for t = dt:dt:T
    u = A\u0(2:end-1);
    u = [0;u;0];
    F = u.*[u(2:end);0] - [0;u(1:end-1)].*[0;u(2:end)];
    u0 = u0 - dt*J*F;
end
figure(1)
plot(x,u0,'linewidth',2)
xlabel('x')
ylabel('u')
title('Legendre-Galerkin')

% Legendre-collocation solver
basis = legendre_basis(N,x);
A = eye(N) - 0.5*dt*mu*L^2;
for t = dt:dt:T
    u = A\(basis'*u0);
    F = u.*[u(2:end);0] - [0;u(1:end-1)].*[0;u(2:end)];
    u0 = basis*J*F;
end
figure(2)
plot(x,u0,'linewidth',2)
xlabel('x')
ylabel('u')
title('Legendre-Collocation')