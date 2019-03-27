function [t, y] = solveEOM(tspan, y0, fs)
% Solve the equation of motion numerically using ode45.
% m (mass) assumed to be 1

N = length(fs.x);

function dydt = odefunc(t, y)
    Phi = reshape(y(1:N^2), N, N);
    dPhi2_dt2 = -FSLaplacian.do(fs, Phi) - Phi;
    dydt = [y(N^2+1:2*N^2); dPhi2_dt2(:)];
end

[t, y] = ode45(@odefunc, tspan, y0);
end