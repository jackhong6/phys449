N = 20;
fs = FuzzySphere(1, N, true);
ss = StringState(pi/2, fs);
M = ss.getM;
iM = ss.getiM;

lat = linspace(-90, 90, 100);
long = linspace(-180, 180, 200);

Z = zeros(length(lat), length(long));

for ii = 1:length(lat)
    theta = deg2rad(lat(ii));
    for jj = 1:length(long)
        phi = deg2rad(long(jj));
        cs = CoherentState([phi, theta, 1], fs, CoordType.spherical);
        z1 = abs(cs.v(:)' * (0.5*(M - 1i*iM)) * cs.v(:));
        %z2 = abs(cs.v(:)' * iM * cs.v(:));
        Z(ii, jj) = z1;
    end
end

[longM, latM] = meshgrid(long, lat);

subplot(2, 1, 1)
axesm('mollweid')
geoshow(latM, longM, Z, 'DisplayType', 'texturemap');
title('Overlap of string state with coherent states ($\langle x | j \rangle \langle -j  | x \rangle$)', 'interpreter', 'latex', 'fontsize', 20)
%caxis([0, 1]);
colorbar;

subplot(2, 1, 2)
ss.draw
title('Overlap of string state with other string states (tr($|a \rangle \langle b |j \rangle \langle -j |$))', 'interpreter', 'latex', 'fontsize', 20)