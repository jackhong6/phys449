[lat, long] = meshgrid(linspace(-90, 90, 20), linspace(-180, 180, 40));
Z = randn(size(lat));
axesm eckert4;
framem; gridm;
h = geoshow(lat, long, Z, 'DisplayType', 'texturemap');
