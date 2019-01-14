classdef Laplacian
    %LAPLACIAN Laplacian on a fuzzy sphere: [Ji, [Ji, Phi]]
    %   Contains the N^2 x N^2 matrix representing the action of minus the
    %   laplacian on the parametrization vector from the StringState class,
    %   as well as the eigenvalues and eigenvectors.
    
    properties
        % K is the N^2 by N^2 matrix representing the action of minus the
        % laplacian on the parametrization vector p (see StringState.m).
        % This matrix is symmetric and has eigenvalues l(l+1)/R^2 for
        % l=0,1,2,3,...,N-1 with multiplicities 2l+1.
        %
        % V is the matrix whose columns are the eigenvectors of K
        % D is the diagonal matrix of eigenvalues corresponding to V
        % R is the radius of the sphere.
        K, V, D, R
    end
    
    methods
       function obj = Laplacian(fs)
            %LAPLACIAN Construct an instance of this class
            %   fs is an instance of FuzzySphere
            %   Initialize K and find V and D
            obj.R = fs.R;
            N = size(fs.x, 1);

            K = zeros(N^2, N^2);
            
            vr = zeros(N^2, 1);
            vc = zeros(N^2, 1);
            
            for r = 1:length(vr)
                vr(r) = 1;
                for c = 1:length(vc)
                    vc(c) = 1;
                    M = StringState.p2M(vc);
                    L = -Laplacian.do(fs, M);
                    p = StringState.M2p(L);
                    K(r, c) = vr' * p;
                    vc(c) = 0;
                end
                vr(r) = 0;
            end
            
            obj.K = sparse(K);
       end
       
% ============= Old code using symbolic programming. Very Slow. ===========
%         function obj = Laplacian(fs)
%             %LAPLACIAN Construct an instance of this class
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
    end
end

