classdef FSLaplacian < handle
    %FSLAPLACIAN Laplacian on a fuzzy sphere: [Li, [Li, Phi]]
    %   Contains the N^2 x N^2 matrix representing the action of minus the
    %   laplacian on the parametrization vector from the StringState class,
    %   as well as the eigenvalues and eigenvectors.
    
    properties
        % K is the cell array of sub-blocks of the N^2 by N^2 matrix 
        % representing the action of minus the laplacian on the 
        % parametrization vector p (see StringState.m).
        % The matrix is bisymmetric and has eigenvalues l(l+1)/R^2 for
        % l=0,1,2,3,...,N-1 with multiplicities 2l+1.
        % There are N bisymmetric sub-blocks of sizes N, 2(N-1), 2(N-2), ...
        % 2(N-(N-1)). 
        %
        % V is the cell array of eigenvectors of the sub-blocks of K
        % D is the cell array of eigenvalues corresponding to V
        % fs is a reference to a FuzzySphere object.
        K, V, D, fs
    end
    
    methods
        function obj = FSLaplacian(fs)
            %LAPLACIAN Construct an instance of this class
            %   fs is an instance of FuzzySphere
            %   Initialize K and find V and D. Use bisymmetric property of K.
            obj.fs = fs;
            
            N = size(fs.x, 1);
            
            obj.K = cell(N, 1);
            obj.V = cell(N, 1);
            obj.D = cell(N, 1);
            
            % Calculate the first block
            obj.K{1} = zeros(N);
            
            % Calculate main diagonal of the first block.
            for ii = 1 : ceil(N/2)
                Kii = obj.calculateKij(ii, ii, N);
                obj.K{1}(ii, ii) = Kii;
                obj.K{1}(N+1-ii, N+1-ii) = Kii;
            end
            
            % Calculate off diagonal (+1/-1) elements of block 1.
            for ii = 1 : ceil((N-1) / 2)
                jj = ii + 1;
                Kij =  obj.calculateKij(ii, jj, N);
                obj.K{1}(ii, jj) = Kij;
                obj.K{1}(jj, ii) = Kij;
                obj.K{1}(N-jj+1, N-ii+1) = Kij;
                obj.K{1}(N-ii+1, N-jj+1) = Kij;
            end
            
            [obj.V{1}, obj.D{1}] = eigs(obj.K{1}, N);
            
            % blocks 2 to N
            offset = N;
            for blk = 2:N
               blk_size = 2 * (N - blk + 1);
               obj.K{blk} = zeros(blk_size); 
               
               % Calculate main diagonal
               for ii = 1 : blk_size / 2
                   Kii = obj.calculateKij(ii+offset, ii+offset, N);
                   obj.K{blk}(ii, ii) = Kii;
                   obj.K{blk}(blk_size+1-ii, blk_size+1-ii) = Kii;
               end
               
               % Calculate the two off diagonals
               for ii = 1 : (blk_size-2) / 2
                   jj = ii + 2;
                   Kij = obj.calculateKij(ii+offset, jj+offset, N);
                   obj.K{blk}(ii, jj) = Kij;
                   obj.K{blk}(jj, ii) = Kij;
                   obj.K{blk}(blk_size-jj+1, blk_size-ii+1) = Kij;
                   obj.K{blk}(blk_size-ii+1, blk_size-jj+1) = Kij;
               end
               
               [obj.V{blk}, obj.D{blk}] = eigs(obj.K{blk}, blk_size);
               offset = offset + blk_size;
            end
            obj.D
        end
        
        function Kij = calculateKij(self, ii, jj, N)
            % Calculate the (i,j) element of the matrix K.
            % v is assumed to be a zero vector.
            v = zeros(N^2, 1);
            v(jj) = 1;
            M = StringState.p2M(v);
            M = -FSLaplacian.do(self.fs, M);
            p = StringState.M2p(M);
            Kij = p(ii);
        end
        
        function Kv = Ktimes(self, v)
            % Return Kv, the result of left multiplying v by K
            N = size(self.fs.x, 1);
            assert(isvector(v) && length(v) == N^2)
            v = v(:);
            Kv = zeros(N^2, 1);
            Kv(1:N) = self.K{1} * v(1:N);
            
            offset = N;
            for blk = 2:N
                blk_size = 2 * (N - blk + 1);
                start_ind = offset+1;
                end_ind = offset + blk_size;
                Kv(start_ind : end_ind) = self.K{blk} * v(start_ind:end_ind);
                offset = offset + blk_size;
            end
        end
        
        function K = getFullK(self)
            K = blkdiag(self.K{:});
        end
        
        function V = getFullV(self)
            V = blkdiag(self.V{:});
        end
        
        function D = getFullD(self)
            D = blkdiag(self.D{:});
        end
       
% ==== Faster code that is still pretty slow. Calculates every element. ===      
%        function obj = Laplacian(fs)
%             %FSLAPLACIAN Construct an instance of this class
%             %   fs is an instance of FuzzySphere
%             %   Initialize K and find V and D
%             obj.R = fs.R;
%             N = size(fs.x, 1);
% 
%             K = zeros(N^2, N^2);
%             
%             vr = zeros(N^2, 1);
%             vc = zeros(N^2, 1);
%             
%             for r = 1:length(vr)
%                 vr(r) = 1;
%                 for c = 1:length(vc)
%                     vc(c) = 1;
%                     M = StringState.p2M(vc);
%                     L = -Laplacian.do(fs, M);
%                     p = StringState.M2p(L);
%                     K(r, c) = vr' * p;
%                     vc(c) = 0;
%                 end
%                 vr(r) = 0;
%             end
%             
%             obj.K = sparse(K);
%        end
% =========================================================================
       
% ============= Old code using symbolic programming. Very Slow. ===========
%         function obj = Laplacian(fs)
%             %FSLAPLACIAN Construct an instance of this class
%             %   fs is an instance of FuzzySphere
%             %   Initialize K and find V and D
%             obj.R = fs.R;
%             N = size(fs.x, 1);
%             x = sym('x', [N^2, 1]);
%             M = StringState.p2M(x);
%             L = simplify(-Laplacian.do(fs, M));
%             p = StringState.M2p(L);
%             obj.K = zeros(N^2, N^2);
%             
%             v = zeros(N^2, 1);
%             for r = 1:length(p)
%                 for c = 1:N^2
%                     v(c) = 1;
%                     obj.K(r, c) = subs(p(r), x, v);
%                     v(c) = 0;
%                 end
%             end
%         end
% =========================================================================
    end
        
    methods (Static)
        function result = do(fs, Phi)
        % LAPLACIAN Return the laplacian of a string state by applying the
        % commutator with the angular momentum matrices [Li, [Li,Phi]]
        N = size(fs.x, 1);
    
        A = - (N^2 - 1) / (4 * fs.R^4);

        result = A * (commutator(fs.x, commutator(fs.x, Phi))...
            + commutator(fs.y, commutator(fs.y, Phi))...
            + commutator(fs.z, commutator(fs.z, Phi)));
        end
               
        function k = p2kBasis(la, p)
            % Change the basis representation of the vector p;
            % la is a Laplacian object.
            % p is in the Lz basis while k is in the basis of the
            % eigenvectors of K.
            N = size(la.fs.x, 1);
            assert(isvector(p) && length(p) == N^2)
            p = p(:);
            k = zeros(N^2, 1);
            k(1:N) = la.V{1}' * p(1:N);
            
            offset = N;
            for blk = 2:N
                blk_size =  2 * (N - blk + 1);
                k(offset+1 : offset+blk_size) = ...
                    la.V{blk}' * p(offset+1 : offset+blk_size);
                offset = offset + blk_size;
            end
        end
        
        function p = k2pBasis(la, k)
            % Change the basis representation of the vector k;
            % p is in the Lz basis while k is in the basis of the
            % eigenvectors of K.
            N = size(la.fs.x, 1);
            assert(isvector(k) && length(k) == N^2)
            k = k(:);
            p = zeros(N^2, 1);
            p(1:N) = la.V{1} * k(1:N);
            
            offset = N;
            for blk = 2:N
                blk_size =  2 * (N - blk + 1);
                p(offset+1 : offset+blk_size) = ...
                    la.V{blk} * k(offset+1 : offset+blk_size);
                offset = offset + blk_size;
            end
        end
    end
end

