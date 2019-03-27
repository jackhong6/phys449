function [F, V] = animateStringState2(t, ss)
%ANIMATESTRINGSTATE Create video of the time evolution of a string state.

N = size(ss.fs.x, 1);
n_theta = 10;
n_phi = 4 * n_theta;

lat = linspace(-90, 90, 2*n_theta);
long = linspace(-180, 180, n_phi);

Z = zeros(size(lat));

fig = figure();
ax = axesm('mollweid');
F(length(t)) = struct('cdata',[],'colormap',[]);

for ii = 1:n_theta
    for jj = 1:n_phi
        ssb = StringState(deg2rad(lat(ii)), deg2rad(long(jj)), ss.fs);
        overlap = StringState.overlap(ss, ssb);
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
    frame = updateFrame(Z, ss, lat, long, latM, longM, n_theta, n_phi, t(n));
    F(n) = frame;
    writeVideo(V, frame);
end

close(V);
end

function frame = updateFrame(Z, ss, lat, long, latM, longM, n_theta, n_phi, t)
for ii = 1:n_theta
    for jj = 1:n_phi
        ssb = StringState(deg2rad(lat(ii)), deg2rad(long(jj)), ss.fs);
        overlap = StringState.overlap(ss.kt(t), ssb.k0, ssb.ik0);
        Z(jj, ii) = overlap;
        Z(jj, 2*n_theta - ii + 1) = overlap;
    end
end

geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
frame = getframe(gcf);

end