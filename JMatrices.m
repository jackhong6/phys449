classdef JMatrices
    % Class containing the three NxN angular momentum operators in the
    % |jm> basis where Jz is diagonal.
    %    Matrices are initialized on object creation.
    %    For example, to create 3x3 matrices:
    %
    %        J = JMatrices(3);
    %
    %    To access the three matrices, you can use
    %
    %        J.x, J.y, J.zx
    %
    %    or equivalently,
    %
    %        J(1), J(2), J(3)
    
    properties
        % x, y, z are NxN angular momentum matrices
        x, y, z
    end
    
    methods
        function J = JMatrices(N)
            % Constructor. Initialize the 3 NxN J.x, J.y, and J.z matrices
            % in the basis where J.z is diagonal.
            j = (N-1) / 2;
            m = j : -1 : -j;
            
            cp = sqrt((j - m(2:end)) .* (j + m(2:end) + 1));
            Jplus = diag(cp, 1);
            
            cm = sqrt((j + m(1:end-1)) .* (j - m(1:end-1) + 1));
            Jminus = diag(cm, -1);
            
            J.x = (Jplus + Jminus) / 2;
            J.y = (Jplus - Jminus) / 2i;
            J.z = diag(m);
        end
        
        function sref = subsref(obj, s)
            % Overload indexing so that J(1) == J.x, J(2) == J.y and
            % J(3) == J.z.
            switch s(1).type
                case '.'
                    sref = builtin('subsref', obj, s);
                case '()'
                    if length(s) == 1
                        if length(s(1).subs{1}) == 1
                            switch(s(1).subs{1})
                                case 1
                                    sref = obj.x;
                                case 2
                                    sref = obj.y;
                                case 3
                                    sref = obj.z;
                                otherwise
                                    error('JMatrices:IndexOutOfRange', 'Index out of range')
                            end
                        else
                            error('JMatrices:IndexingByVectorsNotSupported', 'Index with scalars only')
                        end
                    else
                        % Use built-in for any other expression
                        [sref{1:nargout}] = builtin('subsref',obj,s);
                    end
                case '{}'
                    error('JMatrices:subsref',...
                        'Not a supported subscripted reference')
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        
        function Jn = dot(J, n, coordType)
            % Return Jn, the NxN matrix resulting from the dot product between
            % a unit vector n and the 3 J matrices.
            %
            % ASSUME: n is a 3x1 or 1x3 unit vector. If norm(n) is not 1,
            % then normalize n before applying dot product.
            % If coordType==CoordType.cartesian, n=[x, y, z]
            % If coordType==CoordType.spherical, n=[azimuth, elevation, r]
            %
            % See https://www.mathworks.com/help/matlab/ref/sph2cart.html
            
            switch coordType
                case CoordType.cartesian
                    if norm(n) ~= 1
                        n = n / norm(n);
                    end
                case CoordType.spherical
                    if n(3) ~= 1
                        n(3) = 1;
                    end
                    [n(1), n(2), n(3)] = sph2cart(n(1), n(2), n(3));
            end
            
            Jn = n(1)*J.x + n(2)*J.y + n(3)*J.z;
        end
    end
end