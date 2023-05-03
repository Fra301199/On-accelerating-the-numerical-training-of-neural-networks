function [u_fine, t_fine] = fwd_Euler_parareal(t0, tN, y0, f, dt, dT, DT, k_max)

%[t,y_fine] = fwd_Euler(t0,tN,y0,dt,f);

%% Step 1: coarse approximation for U_0, initial prediction

%DT = 1; % subintervals' length

u0 = y0; % initial condition

%dT = 1;

t_coarse = [t0:dT:tN];
L_coarse = length(t_coarse);

[t_coarse,U_0] = fwd_Euler(t0,tN,u0,dT,f);

% if I don't write explicitely (1:end) the legend is not correct
plot(t_coarse(1:end), U_0(1:end), 'x', 'LineWidth', 2, 'Color', 'r')
legend('U_0')


%% Step 2: iterative steps of parareal

U = U_0;

for m = 1 : (tN - t0) / dT
    t_fine(m,:) = [t0 + (m-1) * dT : dt : t0 + m * dT];
end

L_fine = length(t_fine);

u_fine = zeros(L_coarse, L_fine);
u_fine(1) = u0;

epsilon = 1; % still to implement, for convergence condition
update_rate = zeros(L_coarse, 1);

% implement: min(update_rate(:) < epsilon) == 1;

for k = 1 : k_max

    % Parallel step: fine approximation of the solution: F(tn,tn-1,(U^k)_n-1)
    parfor n = 1 : L_coarse - 1
        [t_fine(n,:),u_fine(n,:)] = fwd_Euler(t_coarse(n),t_coarse(n+1),U_0(n),dt,f);
    end

    figure
    for n = 1:(tN - t0)/dT
        plot(t_fine(n,1:end),u_fine(n,1:end),'-','Linewidth',2, 'Color', 'k')
        hold on
    end

    % Update of the coarse solution
    for h = 1 : L_coarse - 1
        du = sin(t_coarse(h))*U(h) + t_coarse(h);
        U_k(h+1) = U(h) + dT*du;
    end

    % Correction step
    for n = 1 : L_coarse - 1
        du = sin(t_coarse(n))*U(n) + t_coarse(n);
        U(n+1) = U(n) + dT*du;

        U(n+1) = U(n+1) + u_fine(n, end) - U_k(n+1);
    end

    plot(t_coarse, U, 'x', 'LineWidth', 2, 'Color', 'r')
    plot(t_coarse, U_0, 'o', 'LineWidth', 2, 'Color', 'b')
    
    title(['Iteration k = ', num2str(k)]);
    %for i = 1:L_coarse
    %    update_rate(i) = U(i) - U_0(i);
    %end

    U_0 = U;

end

end