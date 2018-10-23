function [t, y] = solveEOM(tspan, y0, m, J)
% Solve the equation of motion numerically using ode45.

N = length(J.x);

    function dydt = odefunc(t, y)
        Phi = reshape(y(1:N^2), [N, N]);
        dPhi_dt = -(laplacian(J, Phi) + m^2 * Phi);
        dydt = [y(N^2+1:2*N^2); dPhi_dt(:)];
    end

[t, y] = ode45(@odefunc, tspan, y0);
end