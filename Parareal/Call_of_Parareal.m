clear all
close all
clc

f = @(t,y) sin(t).*y + t; 
t0 = 0; tN = 14;
y0 = 0;
dt = 0.02;
dT = 1;
DT = dT;
k_max = 8;

[u,t] = fwd_Euler_parareal(t0, tN, y0, f, dt, dT, DT, k_max);

[t_fwd, u_fwd] = fwd_Euler(t0,tN,y0,dt,f);

% Plot of the last iteration of Parareal
% close all
figure
plot(t_fwd,u_fwd,'-','Linewidth',2, 'Color', 'r')
hold on
for n = 1:(tN - t0)/dT
    plot(t(n,1),u(n,1),'o','Linewidth',2, 'Color', 'b')
    plot(t(n,2:end),u(n,2:end),'--','Linewidth',2, 'Color', 'b')
end
legend('Euler','Parareal')