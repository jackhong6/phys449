function [F, V] = animateStringState(t, ss)
%ANIMATESTRINGSTATE Create video of the time evolution of a string state.

N = size(ss.fs.x, 1);
n_theta = 10;
n_phi = 4 * n_theta;

lat = linspace(-90, 90, 2*n_theta);
long = linspace(-180, 180, n_phi);

k01 = zeros(n_theta, n_phi, N^2);
k02 = zeros(n_theta, n_phi, N^2);

for ii = 1:n_theta
    theta = deg2rad(lat(ii));
    for jj = 1:n_phi
        phi = deg2rad(long(jj));
        ssb = StringState(theta, phi, ss.fs);
        k01(ii, jj, :) = ssb.k0;
        k02(ii, jj, :) = ssb.ik0;
    end
end

Z = zeros(size(lat));

fig = figure();
ax = axesm('mollweid');
F(length(t)) = struct('cdata',[],'colormap',[]);

for ii = 1:n_theta
    for jj = 1:n_phi
        overlap1 = 2*abs(dot(squeeze(k01(ii, jj, :)), ss.kt(t(1))));
        overlap2 = 2*abs(dot(squeeze(k02(ii, jj, :)), ss.kt(t(1))));
        overlap = sqrt(overlap1^2 + overlap2^2);
        Z(jj, ii) = overlap;
        Z(jj, 2*n_theta - ii + 1) = overlap;
    end
end

[latM, longM] = meshgrid(lat, long);

geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
caxis([0 1]);
colorbar;
title('Time evolution of StringState')

F(1) = getframe(gcf);
V = VideoWriter(['../Videos/string_state_animation_N' num2str(N) '.mp4'], 'MPEG-4');
V.FrameRate = 7;
open(V);

for n = 2:length(t)
    kt = ss.kt(t(n));
    frame = updateFrame(Z, k01, k02, kt, latM, longM, n_theta, n_phi);
    F(n) = frame;
    writeVideo(V, frame);
end

close(V);
end

function frame = updateFrame(Z, k01, k02, kt, latM, longM, n_theta, n_phi)
for ii = 1:n_theta
    for jj = 1:n_phi
        overlap1 = 2*abs(dot(squeeze(k01(ii, jj, :)), kt));
        overlap2 = 2*abs(dot(squeeze(k02(ii, jj, :)), kt));
        overlap = sqrt(overlap1^2 + overlap2^2);
        Z(jj, ii) = overlap;
        Z(jj, 2*n_theta - ii + 1) = overlap;
    end
end

geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
frame = getframe(gcf);

end