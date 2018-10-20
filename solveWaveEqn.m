function [t, y] = solveWaveEqn(t, m, J, Phi, dPhi_dt0)
% Solve the wave equation numerically using ode45.

    function dydt = odefunc(t, y)
        dydt = [y(2); J.laplacian(];
    end
end



