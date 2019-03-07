classdef CoherentState
    % COHERENTSTATE Nx1 vector representation of a coherent state in Lz basis.
    %   A coherent state is the eigenvector with the largest eigenvalue
    %   of the angular momentum operator in the n direction.
    
    properties
        % v is the Nx1 eigenvector of Ln with eigenvalue j in the Lz basis, 
        % where Ln is the angular momentum operator in the n direction.
        v, n, coordType
    end
    
    methods
        function obj = CoherentState(n, fs, coordType)
            % CONSTRUCTOR
            % If coordType==CoordType.cartesian, n=[x, y, z]
            % If coordType==CoordType.spherical, n=[azimuth, elevation, r]
            if nargin > 0
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
    end
end