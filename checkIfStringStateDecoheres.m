N = 10;
fs = FuzzySphere(1, N, true);

ss = StringState(pi/4, fs);

t = linspace(0, 10, 500);
x = zeros(1, length(t));

for k = 1:length(t)
    Mt1 = ss.getMt(t(k));
    Mt2 = ss.getiMt(t(k));
    A = 0.5 * (Mt1 - 1i*Mt2);
    x(k) = trace((A'*A)^2 - A'*A)^2;
end

plot(t, x)

dPhi_dt_0 = zeros(N);
y0 = [Phi.M(:); dPhi_dt_0(:)];
[t, y] = solveEOM(tspan, y0, m, J);