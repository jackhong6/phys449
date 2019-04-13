function [F, V] = animateStringState2(t, ss)
%ANIMATESTRINGSTATE2 Create video of the time evolution of a string state 
% using horizontally symmetric states.

N = size(ss.fs.x, 1);
n_phi = 20;
n_theta = n_phi;

lat = linspace(-90, 90, n_theta);
long = linspace(-180, 180, 2*n_phi);

k01 = zeros(n_phi, n_theta, N^2);
k02 = zeros(n_phi, n_theta, N^2);

for ii = 1:n_theta
    theta = deg2rad(lat(ii));
    for jj = 1:n_phi
        phi = deg2rad(long(jj));
        n1 = [phi, theta, 1];
        n2 = [-phi, theta, 1];
        cs1 = CoherentState(n1, ss.fs);
        cs2 = CoherentState(n2, ss.fs);
        ssb = StringState(cs1, cs2, ss.fs);
        k01(jj, ii, :) = ssb.k0;
        k02(jj, ii, :) = ssb.ik0;
    end
end

[longM, latM] = meshgrid(long, lat);
Z = zeros(size(latM));

figure();
axesm('mollweid');
F(length(t)) = struct('cdata',[],'colormap',[]);

for ii = 1:n_theta
    for jj = 1:n_phi
        overlap1 = 2*abs(dot(squeeze(k01(jj, ii, :)), ss.kt(t(1))));
        overlap2 = 2*abs(dot(squeeze(k02(jj, ii, :)), ss.kt(t(1))));
        overlap = sqrt(overlap1^2 + overlap2^2);
        Z(ii, jj) = overlap;
        Z(ii, 2*n_phi - jj + 1) = overlap;
    end
end

geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
caxis([0 1]);
colorbar;
title('Time evolution of StringState')

F(1) = getframe(gcf);
V = VideoWriter(['../Videos/animateStringState2_N' num2str(N) '.mp4'], 'MPEG-4');
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
        overlap1 = 2*abs(dot(squeeze(k01(jj, ii, :)), kt));
        overlap2 = 2*abs(dot(squeeze(k02(jj, ii, :)), kt));
        overlap = sqrt(overlap1^2 + overlap2^2);
        Z(ii, jj) = overlap;
        Z(ii, 2*n_phi - jj + 1) = overlap;
    end
end

geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
frame = getframe(gcf);

end