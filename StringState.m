classdef StringState
    % STRINGSTATE A string state is the outer product of two coherent states.
    %   It is usually represented by an NxN hermitian matrix, but is stored
    %   as a parametrization vector p described below.
    
    properties
        % p is the parametrization of the string state by a vector of length N^2.
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
        % fs is the FuzzySphere that the string state lives on
        % m is the mass of the StringState
        p, fs, m
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
                    cs1 = CoherentState(obj.fs, n1, CoordType.spherical);
                    cs2 = CoherentState(obj.fs, n2, CoordType.spherical);
                    obj = StringState(cs1, cs2, obj.fs);
                else
                    % Interpret as p vector 
                    obj.p = varargin{1};
                end
            elseif nargin == 3
                % If there are 3 arguments, then assume the first 2 are
                % of type CoherentState.
                obj.fs = varargin{3};

                M = (varargin{1}.v(:) * varargin{2}.v(:)') ...
                    + (varargin{2}.v(:) * varargin{1}.v(:)');
                
                % Normalize the string state so that Tr(A' * A) = 1
                A = trace(M' * M);
                M = M / sqrt(A);
                
                obj.p = obj.M2p(M);        
            else
                error('Inputs: ([CoherentStates, parametrization vector, opening angle]; FuzzySphere)')
            end
            obj.m = 1;
        end
        
        function M = getM(self)
            M = self.p2M(self.p);
        end
        
        function k = getk(self)
            k = FSLaplacian.p2kBasis(self.fs.la, self.p);
        end
        
        function k_t = kt(self, t, k0, dkdt0)
            w = sqrt(diag(self.fs.la.getFullD) + self.m^2);
            k_t = k0.*cos(w*t) + (dkdt0 ./ w) .* sin(w*t);
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