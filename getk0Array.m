function k0 = getk0Array(fs, inputArg2)
%GETk0ARRAY Summary of this function goes here
%   Detailed explanation goes here
N = size(fs.x, 1);
n_theta = 50;
n_phi = 4 * n_theta;

[lat, long] = meshgrid(linspace(-90, 90, 2*n_theta), linspace(-180, 180, n_phi));

k0 = zeros(n_theta, n_phi, N^2);
for ii = 1:n_theta
    theta = deg2rad(lat(1, ii));
    for jj = 1:n_phi
        phi = deg2rad(long(jj, 1));
        ss = StringState(theta, phi, fs);
        k0(ii, jj, :) = ss.calculate_k0;
    end
end
end

