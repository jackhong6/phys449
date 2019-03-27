N = 100;
fs = FuzzySphere(1, N, true);

theta = deg2rad(45);
ss = StringState(theta, fs);

t = linspace(0, 20, 200);
x = zeros(1, length(t));
%Aeigs = zeros(length(t), N);

% M0 = ss.getM;
% iM0 = ss.getiM;
% A = 0.5 * (M0 - 1i*iM0);
% C = trace(A' * A);

for k = 1:length(t)
    Mt1 = ss.getMt(t(k));
    Mt2 = ss.getiMt(t(k));
    A = 0.5 * (Mt1 - 1i*Mt2);
    %Aeigs(k, :) = eigs(A, N);
    C = trace(A' * A);
    x(k) = trace((A'*A)^2 / C^2 - (A'*A) / C)^2;
end
figure()
plot(t, x)
title(['N = ' num2str(N) ' ; opening angle = ' num2str(rad2deg(theta))])
xlabel('Time', 'interpreter', 'latex')
ylabel('tr$((A^\dagger A)^2 - A^\dagger A)^2$', 'interpreter', 'latex')

% tspan = [0, 3];
% v0 = zeros(N);
% M1 = ss.getM;
% M2 = ss.getiM; 
% y01 = [M1(:); v0(:)];
% y02 = [M2(:); v0(:)];
% [t1, y1] = solveEOM(tspan, y01, fs);
% [t2, y2] = solveEOM(tspan, y02, fs);
% x = zeros(1, length(t1));
% 
% for k = 1:length(t1)
%     Mt1 = reshape(y1(k, 1:N^2), N, N);
%     Mt2 = reshape(y2(k, 1:N^2), N, N);
%     A = 0.5 * (Mt1 - 1i*Mt2);
%     x(k) = trace((A'*A)^2 - A'*A)^2;
% end
% 
% figure()
% plot(t1, x);
% ylim([0, 0.07]);