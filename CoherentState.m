classdef CoherentState
    % COHERENTSTATE Nx1 vector representation of a coherent state in Jz basis.
    %   A coherent state is the eigenvector with the largest eigenvalue
    %   of the angular momentum operator in the n direction.
    
    properties
        % v is the Nx1 eigenvector of Ln with eigenvalue j, 
        % where Ln is the angular momentum operator in the n direction.
        v, n, coordType
    end
    
    methods
        function obj = CoherentState(fs, n, coordType)
            % CONSTRUCTOR
            % If coordType==CoordType.cartesian, n=[x, y, z]
            % If coordType==CoordType.spherical, n=[azimuth, elevation, r]
            if nargin > 0
                N = size(fs.x, 1);
                obj.n = n;
                obj.coordType = coordType;
                j = (N - 1) / 2;
                eig_value = sqrt(j / (j + 1)) * fs.R;

                tol = 1e-9;
                Xn = fs.dot(n, coordType);

                [V, D] = eigs(Xn, 2);
                obj.v = V(:, abs(diag(D) - eig_value) < tol);
            end
        end
    end
end