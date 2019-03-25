addpath('subaxis')

figure()

fs = FuzzySphere(1, 5, false);
cs = CoherentState([0, 0, 1], fs, CoordType.spherical);
subaxis(2, 2, 1, 'm', 0, 'p', 0, 's', 0, 'holdaxis', 1);
%ax1.OuterPosition = [0, 0.5, 0.5, 1];
cs.draw
colorbar off

fs = FuzzySphere(1, 10, false);
cs = CoherentState([0, 0, 1], fs, CoordType.spherical);
subaxis(2, 2, 2, 'm', 0, 'p', 0, 's', 0, 'mr', 0.1, 'holdaxis', 1);
%ax2.OuterPosition = [0.5, 0.5, 0.5, 1];
cs.draw
colorbar off

fs = FuzzySphere(1, 20, false);
cs = CoherentState([0, 0, 1], fs, CoordType.spherical);
subaxis(2, 2, 3, 'm', 0, 'p', 0, 's', 0, 'holdaxis', 1);
%ax3.OuterPosition = [0, 0, 0.5, 1];
cs.draw
colorbar off

fs = FuzzySphere(1, 40, false);
cs = CoherentState([0, 0, 1], fs, CoordType.spherical);
subaxis(2, 2, 4, 'm', 0, 'p', 0, 's', 0, 'mr', 0.1, 'holdaxis', 1);
%ax4.OuterPosition = [0.5, 0, 0.5, 1];
cs.draw
colorbar off
colorbar('manual', 'position', [0.9, 0.05, .03, 0.9])
