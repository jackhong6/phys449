classdef StringState
    % STRINGSTATE A string state is the outer product of two coherent states.
    
    properties (Access = private)
        % M is the matrix representation in the Jz basis.
        % p is the parametrization of M by a vector of length N^2.
        % These properties are private to ensure that they are always in
        % sync and represent the same string state.
        %
        % Example: for a 3x3 string state the parametrization is
        %      p                          M
        %     _  _
        %    | x1 |
        %    | x2 |       _                                  _
        %    | x3 |      | sqrt(2) x1, x4 + ix5  , x8 + ix9   |
        %    | x4 |      |                                    |
        %    | x5 | <--> | x4 - ix5  , sqrt(2) x2, x6 + ix7   |
        %    | x6 |      |                                    |
        %    | x7 |      | x8 - ix9  , x6 - ix7  , sqrt(2) x3 |
        %    | x8 |       ?                                  ?
        %    | x9 |
        %     ?  ?
        M, p
    end
    
    methods
        function obj = StringState(varargin)
            % CONSTRUCTOR Create an instance of this class using either two
            % CoherentStates or a p vector parametrization.
            if nargin == 2
                % If there are two arguments, then assume they are of type 
                % CoherentState 
                obj.M = (varargin{1}.v * varargin{2}.v') ...
                    + (varargin{2}.v * varargin{1}.v');
                
                % Normalize the string state so that Tr(A' * A) = 1
                A = trace(obj.M' * obj.M);
                obj.M = obj.M / sqrt(A);
                
                obj.p = obj.M2p(obj.M);
                
            elseif nargin == 1
                obj.M = p2M(varargin{1});
                obj.p = varargin{1};
            else
                return
            end
        end
        
        function M = getM(self)
            M = self.M;
        end
    end
    
    methods (Static)
        function M = p2M(p)
            % P2M Return the matrix corresponding to the parametrization p.
            N = sqrt(length(p));
            M = diag(sqrt(2) * p(1:N));
            start_ind = N + 1;
            
            for k = 1:N-1
                diag_length = N - k;
                end_ind = start_ind + 2*diag_length - 1;
                real_ind = start_ind : 2 : end_ind - 1;
                imag_ind = start_ind + 1 : 2 : end_ind;
                
                M(diag(true(diag_length, 1),  k)) = ...
                    p(real_ind) + 1i*p(imag_ind);

                M(diag(true(diag_length, 1), -k)) = ...
                    p(real_ind) - 1i*p(imag_ind);
                
                start_ind = end_ind + 1;
            end
            
        end
        
        function p = M2p(M)
            % M2P Return the vector parametrization of the matrix M.
            %assert(ishermitian(M));
            N = size(M, 1);
            p = zeros(N^2, 1);
            p(1:N) = diag(M) / sqrt(2);
            start_ind = N + 1;
            
            for k = 1:N-1
                diag_length = N - k;
                end_ind = start_ind + 2*diag_length - 1;
                
                real_ind = start_ind : 2 : end_ind - 1;
                imag_ind = start_ind + 1 : 2 : end_ind;
                
                kth_diag = M(diag(true(diag_length, 1),  k));
                p(real_ind) = real(kth_diag);
                p(imag_ind) = imag(kth_diag);
                
                start_ind = end_ind + 1;
            end
        end
    end
end