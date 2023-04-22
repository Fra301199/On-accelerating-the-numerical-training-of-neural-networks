% Parareal: first implementation attempt
clear all
close all
clc

% Parareal algorithm to solve the ODE
% u'(t) = sen(t)*u(t) + t , t in [0,14]
% u(0) = 1

%--------------------------------------------------------------------------
% Forward Euler in the fine time-mesh
%--------------------------------------------------------------------------

t0 = 0; tN = 14;
dt = 0.02;
N = (tN - t0) / dt;

t = linspace(t0,tN,N+1);
% times = [t0:dT:TN]

y0 = 1;
y_fine = [];
y_fine(1) = y0;

for k = 1:N
    dy = sin(t(k))*y_fine(k) + t(k);
    y_fine(k+1) = y_fine(k) + dt*dy;
end

%plot(t, y_fine, '-', 'LineWidth', 2)

%--------------------------------------------------------------------------
% PARAREAL: più real che para
%--------------------------------------------------------------------------

% DT ~ suppose to want N subintervals => DT = (tN - t0) / N
%
% u_coarse = solution provided by coarse propagator
% dT = coarse time step
%
% u_fine = solution provided by fine propagation operator
% dt = fine time step
%
% U ~ solution from the "Correction Step"

% Forward Euler
t0 = 0; tN = 14;
DT = 1;
u0 = 1;


% Step 1: coarse approximation for U_0, initial prediction
% e.g. Forward Euler with time step dT
dT = 1;
t_coarse = [t0:dT:tN];
L_coarse = length(t_coarse);

U_0 = [];
U_0(1) = u0;

for k = 1 : L_coarse - 1
    du = sin(t_coarse(k))*U_0(k) + t_coarse(k);
    U_0(k+1) = U_0(k) + dT*du;
end

plot(t_coarse, U_0, 'x', 'LineWidth', 2, 'Color', 'r')
hold on
plot(t, y_fine, '-', 'LineWidth', 2, 'Color', 'g')
legend('U^0', 'Fu^0')

% Step 2: iterative step of parareal

%U = [];
%U(1) = U_0(1);
U = U_0;

dt = 0.02;
for m = 1 : (tN - t0)/dT
    t_fine(m,:) = [t0 + (m-1)*dT:dt:t0+m*dT];
end
L_fine = length(t_fine);
N_fine = 1/dt; % n° of "fine subintervals" of each "coarse subint."

u_fine = [];
u_fine(1) = u0;

% very few iterations are necessary for the convergence
%----------
k_max = 4;
%k_max = 10;
%----------

epsilon = 1; % still to implement, for convergence condition

for k = 1 : k_max
    
    figure()
    plot(t, y_fine, '-', 'LineWidth', 2, 'Color', 'g')
    % Parallel step: fine approximation of the solution

    % F(tn,tn-1,(U^k)_n-1)
    for n = 1 : L_coarse - 1
        % I.C. in each subinterval, given by the coarse approx
        u_fine(n, 1) = U_0(n);

        for h = 1 : N_fine
            du = sin(t_fine(n,h))*u_fine(n,h) + t_fine(n,h);
            u_fine(n, h+1) = u_fine(n, h) + du*dt;
        end
    end

    hold on
    for n = 1:(tN - t0)/dT
        plot(t_fine(n,1:50),u_fine(n,1:50),'-','Linewidth',2, 'Color', 'k')
    end

    % update of the coarse solution
    for k = 1 : L_coarse - 1
        du = sin(t_coarse(k))*U(k) + t_coarse(k);
        U_k(k+1) = U(k) + dT*du;
    end


    for n = 1 : L_coarse - 1
        du = sin(t_coarse(n))*U(n) + t_coarse(n);
        U(n+1) = U(n) + dT*du;
        % Correction step
        U(n+1) = U(n+1) + u_fine(n, end) - U_k(n+1);
    end

    plot(t_coarse, U, 'x', 'LineWidth', 2, 'Color', 'r')
    plot(t_coarse, U_0, 'o', 'LineWidth', 2, 'Color', 'b')

    U_0 = U;

    %for k = 1 : L_coarse - 1
    %    du = sin(t_coarse(k))*U_0(k) + t_coarse(k);
    %    U_0(k+1) = U(k) + dT*du;
    %end

end