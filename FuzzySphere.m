classdef FuzzySphere
    % FUZZYSPHERE The fuzzy sphere of radius R with N quantum cells
    % represented by three NxN angular momentum matrices in the |jm> basis 
    % where Jz is diagonal.
    %
    %    For example, to create the fuzzy sphere of radius 5 with 3 cells:
    %
    %        fs = FuzzySphere(5, 3);
    %
    %    To access the three matrices, you can use
    %
    %        fs.x, fs.y, fs.z
    %
    %    or equivalently,
    %
    %        fs(1), fs(2), fs(3)
    
    properties
        % x, y, z are NxN matrices. R is the radius of the sphere.
        x, y, z, R
    end
    
    methods
        function fs = FuzzySphere(R, N)
            % CONSTRUCTOR Initialize the 3 NxN matrices
            % in the basis where fs.z is diagonal.
            if nargin > 0
                j = (N - 1) / 2;
                m = j : -1 : -j;

                A = R / sqrt(j * (j + 1));

                cp = sqrt((j - m(2:end)) .* (j + m(2:end) + 1));
                Jplus = diag(cp, 1);

                cm = sqrt((j + m(1:end-1)) .* (j - m(1:end-1) + 1));
                Jminus = diag(cm, -1);

                fs.R = R;
                fs.x = A * (Jplus + Jminus) / 2;
                fs.y = A * (Jplus - Jminus) / 2i;
                fs.z = A * diag(m);
            end
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
                                    error('FuzzySphere:IndexOutOfRange', 'Index out of range')
                            end
                        else
                            error('FuzzySphere:IndexingByVectorsNotSupported', 'Index with scalars only')
                        end
                    else
                        % Use built-in for any other expression
                        [sref{1:nargout}] = builtin('subsref',obj,s);
                    end
                case '{}'
                    error('FuzzySphere:subsref',...
                        'Not a supported subscripted reference')
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        
        function Xn = dot(fs, n, coordType)
            % FUZZYSPHERE.DOT Return the NxN matrix resulting from the dot
            % product between a unit vector n and the 3 x,y,z matrices.
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
            
            Xn = n(1)*fs.x + n(2)*fs.y + n(3)*fs.z;
        end
    end
end