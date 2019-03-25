n_angles = 300;
azimuths = linspace(-pi, pi, n_angles);
Ns = [5, 10, 20, 40];
overlaps = zeros(length(Ns), n_angles);

for ii = 1:length(Ns)
    fs = FuzzySphere(1, Ns(ii), true);
    n = [0, 0, 1];
    cs = CoherentState(n, fs, CoordType.spherical);
    for jj = 1:length(azimuths)
        n = [azimuths(jj), 0, 1];
        cs1 = CoherentState(n, fs, CoordType.spherical);
        overlaps(ii, jj) = abs(cs1.v(:)' * cs.v(:));
    end
end

hold on
for ii = 1:length(Ns)
    txt = ['N = ',num2str(Ns(ii))];
    plot(rad2deg(azimuths), overlaps(ii, :), ...
    'linewidth', 1.2, 'DisplayName', txt)
end

title('Overlap with coherent states along the equator')
xlabel('Azimuthal angle (degrees)', 'interpreter', 'latex')
ylabel('Overlap $\langle x | p \rangle$', 'interpreter', 'latex')
legend show
