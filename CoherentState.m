classdef CoherentState
    % COHERENTSTATE Nx1 vector representation of a coherent state in Lz basis.
    %   A coherent state is the eigenvector with the largest eigenvalue
    %   of the angular momentum operator in the n direction.
    
    properties
        % v is the Nx1 eigenvector of Ln with eigenvalue j in the Lz basis, 
        % where Ln is the angular momentum operator in the n direction.
        % fs is the FuzzySphere where the coherent state lives on.
        v, n, coordType, fs
    end
    
    methods
        function obj = CoherentState(n, fs, coordType)
            % CONSTRUCTOR
            % If coordType==CoordType.cartesian, n=[x, y, z]
            % If coordType==CoordType.spherical, n=[azimuth, elevation, r]
            if nargin > 0
                obj.fs = fs;
                N = size(fs.x, 1);
                obj.n = n;
                obj.coordType = coordType;

                Xn = fs.dot(n, coordType);

                [V, ~] = eigs(Xn + fs.R*eye(N), 1);
                %obj.v = V(:, abs(diag(D) - eig_value) < tol);
                if sign(real(V(1))) == -1
                    V = -V;
                end
                obj.v = V;
            end
        end
        
        function h = draw(self, varargin)
            % Draw the CoherentState and return the axes. Optional input of
            % two vectors:
            %   lat - latitudes (in degrees) to plot
            %   long - longitudes (in degrees) to plot
            if length(varargin) == 2
                lat = varargin{1};
                long = varargin{2};
                n_theta = length(lat);
                n_phi = length(long);
            else
                n_theta = 100;
                n_phi = 2 * n_theta;
                lat = linspace(-90, 90, n_theta);
                long = linspace(-180, 180, n_phi);
            end
            
            overlaps = zeros(n_theta, n_phi);
            
            for ii = 1:n_theta
                theta = deg2rad(lat(ii));
                for jj = 1:n_phi
                    phi = deg2rad(long(jj));
                    nb = [phi, theta, 1];
                    csb = CoherentState(nb, self.fs, CoordType.spherical);
                    overlaps(ii, jj) = abs(csb.v(:)' * self.v(:));
                end
            end
            
            axesm('mollweid');
            [longM, latM] = meshgrid(long, lat);
            h = geoshow(latM, longM, overlaps, 'DisplayType', 'texturemap');
            caxis([0, 1])
            colorbar;
            %dim = [.7 .4 .3 .4];
            %str = ['N=' num2str(size(self.fs.x, 1))];
            %annotation('textbox',dim,'String',str,'FitBoxToText','on', 'fontsize', 16)
            %title_txt = ['N = ' num2str(size(self.fs.x, 1))];
            %title(title_txt)
        end
    end
end