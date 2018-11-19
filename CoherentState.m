classdef CoherentState
    % Class that represents a coherent state.
    %   A coherent state is the eigenvector with the largest eigenvalue (j)
    %   of the Jn matrix.
    
    properties
        % v is the Nx1 eigenvector of Jn with eigenvalue j, where Jn is the
        % angular momentum operator in the n direction.
        v, n, coordType
    end
    
    methods
        function obj = CoherentState(J, n, coordType)
            % Constructor:
            N = size(J.x, 1);
            obj.n = n;
            obj.coordType = coordType;
            j = (N - 1) / 2;
            
            tol = 1e-9;
            Jn = J.dot(n, coordType);
            
            [V, D] = eigs(Jn, 2);
            obj.v = V(:, abs(diag(D)-j) < tol);
        end
    end
end