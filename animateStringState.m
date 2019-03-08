function F = animateStringState(fs, opening_angle, t)
%ANIMATESTRINGSTATE Create video of the time evolution of a string state.

N = size(fs.x, 1);
n_theta = 50;
n_phi = 4 * n_theta;

[lat, long] = meshgrid(linspace(-90, 90, 2*n_theta), linspace(-180, 180, n_phi));

k01 = zeros(n_theta, n_phi, N^2);
k02 = zeros(n_theta, n_phi, N^2);
for ii = 1:n_theta
    theta = deg2rad(lat(1, ii));
    for jj = 1:n_phi
        phi = deg2rad(long(jj, 1));
        ss1 = StringState(theta, phi, fs);
        ss2 = StringState(theta, phi, fs, 1);
        k01(ii, jj, :) = ss1.calculate_k0;
        k02(ii, jj, :) = ss2.calculate_k0;
    end
end

ss = StringState(opening_angle, fs);
Z = zeros(size(lat));

F(length(t)) = struct('cdata',[],'colormap',[]);

for n = 1:length(t)
    kt = ss.kt(t(n));
    F(n) = makeFrame(Z, k01, k02, kt, lat, long, n_theta, n_phi);
end

end

function frame = makeFrame(Z, k01, k02, kt, lat, long, n_theta, n_phi)
for ii = 1:n_theta
    for jj = 1:n_phi
        overlap1 = 2*abs(dot(squeeze(k01(ii, jj, :)), kt));
        overlap2 = 2*abs(dot(squeeze(k02(ii, jj, :)), kt));
        overlap = overlap1^2 + overlap2^2;
        Z(jj, ii) = overlap;
        Z(jj, 2*n_theta - ii + 1) = overlap;
    end
end
axesm eckert4;
framem; gridm;
geoshow(lat, long, Z, 'DisplayType', 'texturemap');
colorbar;
title('Time evolution of StringState with 45deg opening angle (N=10, v=0)')
frame = getframe(gcf);
end