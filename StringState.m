classdef StringState
    % STRINGSTATE A string state is the outer product of two coherent states.
    %   It is usually represented by an NxN hermitian matrix, but is stored
    %   as a parametrization vector p described below.
    
    properties
        % p is the parametrization of the string state by a vector of length N^2.
        % fs is the FuzzySphere that the string state lives on
        % k0 is the parametrization vector p in the basis of the laplacian
        %   K at time = 0
        % dkdt0
        % m is the mass of the StringState
        % w is a vector of time evolution frequencies for the p vector
        %   elements in the k basis.
        %
        % Example: for a 3x3 string state the parametrization is
        %
        %      p                          M
        %     _  _
        %    | x1 |
        %    | x2 |       _                                  _
        %    | x3 |      | sqrt(2) x1, x4 + ix5  , x8 + ix9   |
        %    | x4 |      |                                    |
        %    | x5 | <--> | x4 - ix5  , sqrt(2) x2, x6 + ix7   |
        %    | x6 |      |                                    |
        %    | x7 |      | x8 - ix9  , x6 - ix7  , sqrt(2) x3 |
        %    | x8 |       -                                  -
        %    | x9 |
        %     -  -
        p, fs, k0, dkdt0, m, w
    end
    
    methods
        function obj = StringState(varargin)
            % CONSTRUCTOR Create an instance of this class using either two
            % CoherentStates, a p vector parametrization, or an opening angle.
            if nargin == 2
                obj.fs = varargin{2};
                if isscalar(varargin{1})
                    % Interpret as an opening angle in radians
                    n1 = [0, varargin{1}, 1];  % [azimuth, elevation, radius]
                    n2 = [0, -varargin{1}, 1];
                    cs1 = CoherentState(n1, obj.fs, CoordType.spherical);
                    cs2 = CoherentState(n2, obj.fs, CoordType.spherical);
                    obj = StringState(cs1, cs2, obj.fs);
                else
                    % Interpret as p vector
                    obj.p = varargin{1};
                end
                
            elseif nargin == 3
                 obj.fs = varargin{3};
                if isa(varargin{1}, 'CoherentState')
                    % Interpret as two CoherentStates
                    M = (varargin{1}.v(:) * varargin{2}.v(:)') ...
                        + (varargin{2}.v(:) * varargin{1}.v(:)');
                    
                    % Normalize the string state so that Tr(M' * M) = 1
                    A = trace(M' * M);
                    M = M / sqrt(A);
                    obj.p = obj.M2p(M);
                else
                    % Interpret as (opening angle, azimuthal angle, etc.)
                    n1 = [varargin{2}, varargin{1}, 1];  % [azimuth, elevation, radius]
                    n2 = [varargin{2}, -varargin{1}, 1];
                    cs1 = CoherentState(n1, obj.fs, CoordType.spherical);
                    cs2 = CoherentState(n2, obj.fs, CoordType.spherical);
                    obj = StringState(cs1, cs2, obj.fs);
                end
            elseif nargin == 4
                obj.fs = varargin{3};
                if isa(varargin{1}, 'CoherentState')
                    % Interpret as two CoherentStates
                    M = 1i*(varargin{1}.v(:) * varargin{2}.v(:)') ...
                        - 1i*(varargin{2}.v(:) * varargin{1}.v(:)');
                    
                    % Normalize the string state so that Tr(M' * M) = 1
                    A = trace(M' * M);
                    M = M / sqrt(A);
                    obj.p = obj.M2p(M);
                else
                    % Interpret as (opening angle, azimuthal angle, etc.)
                    n1 = [varargin{2}, varargin{1}, 1];  % [azimuth, elevation, radius]
                    n2 = [varargin{2}, -varargin{1}, 1];
                    cs1 = CoherentState(n1, obj.fs, CoordType.spherical);
                    cs2 = CoherentState(n2, obj.fs, CoordType.spherical);
                    obj = StringState(cs1, cs2, obj.fs, true);
                end
            else
                error('Inputs: ([CoherentStates, parametrization vector, opening angle]; FuzzySphere)')
            end
            
            obj.m = 1;
            
            if ~isempty(obj.fs.la)
                obj.w = obj.getw();
                obj.k0 = obj.calculate_k0;
                obj.dkdt0 = zeros(length(obj.p), 1);
            end
        end
        
        function overlap = overlap0(self, ss, t)
            % Return the normalized overlap between itself (at t=0) and a
            % StringState ss at time=t
            kt = ss.kt(t);
            overlap = 2*abs(kt(:)' * self.k0(:));
        end
        
        function M = getM(self)
            M = self.p2M(self.p);
        end
        
        function w = getw(self)
            w = sqrt(diag(self.fs.la.getFullD) + self.m^2);
        end
        
        function k0 = calculate_k0(self)
            k0 = FSLaplacian.p2kBasis(self.fs.la, self.p);
        end
        
        function dkdt = calculate_dkdt0(self, v, n, coordType)
            % Return a vector of velocities for the k vector given
            % initial speed v and axis of rotation n.
            % n is assumed to be a unit vector pointing in the direction
            % of the axis of rotation. The direction of rotation is given
            % by the right hand rule.
            dMdt = 1i * v * commutator(self.fs.dot(n, coordType), self.getM);
            dpdt = StringState.M2p(dMdt);
            dkdt = FSLaplacian.p2kBasis(self.fs.la, dpdt);
        end
        
        function k_t = kt(self, t)
            % Return the time evolution k0 at time t. t must be a scalar.
            % Must assign properties w, k0, and dkdt0 before using.
            k_t = self.k0.*cos(self.w*t) + (self.dkdt0 ./ self.w) .* sin(self.w*t);
        end
    end
    
    methods (Static)
        function M = p2M(p)
            % P2M Return the matrix corresponding to the parametrization p.
            N = sqrt(length(p));
            M = diag(sqrt(2) * p(1:N));
            start_ind = N + 1;
            
            for m = 1:N-1
                diag_length = N - m;
                end_ind = start_ind + 2*diag_length - 1;
                real_ind = start_ind : 2 : end_ind - 1;
                imag_ind = start_ind + 1 : 2 : end_ind;
                
                M(diag(true(diag_length, 1),  m)) = ...
                    p(real_ind) + 1i*p(imag_ind);
                
                M(diag(true(diag_length, 1), -m)) = ...
                    p(real_ind) - 1i*p(imag_ind);
                
                start_ind = end_ind + 1;
            end
        end
        
        function p = M2p(M)
            % M2P Return the vector parametrization of the matrix M.
            % assert(ishermitian(M));
            N = size(M, 1);
            
            p = zeros(N^2, 1);
            
            p(1:N) = diag(M) / sqrt(2);
            start_ind = N + 1;
            
            for m = 1:N-1
                diag_length = N - m;
                end_ind = start_ind + 2*diag_length - 1;
                
                real_ind = start_ind : 2 : end_ind - 1;
                imag_ind = start_ind + 1 : 2 : end_ind;
                
                kth_diag = M(diag(true(diag_length, 1),  m));
                p(real_ind) = real(kth_diag);
                p(imag_ind) = imag(kth_diag);
                
                start_ind = end_ind + 1;
            end
        end
    end
end