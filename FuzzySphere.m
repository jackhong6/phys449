classdef FuzzySphere < handle
    % FUZZYSPHERE The fuzzy sphere of radius R with N quantum cells
    % represented by three NxN angular momentum matrices in the |jm> basis 
    % where Lz is diagonal.
    %
    %    For example, to create the fuzzy sphere of radius 2 with 3 cells:
    %
    %        fs = FuzzySphere(2, 3);
    %
    %    To access the three matrices, you can use
    %
    %        fs.x, fs.y, fs.z
    %
    %    or equivalently,
    %
    %        fs(1), fs(2), fs(3)
    
    properties
        % x, y, z are NxN matrices.
        % R is the radius of the sphere.
        % la is the FSLaplacian on the sphere.
        % kv is an array of vertical StringState k vectors on the sphere.
        %     i.e. same longitudes, +/- latitudes. M = |a><b| + |b><a|
        % ikv is the corresponding antisymmetric vectors of kv
        %     iM = i|a><b| - i |b><a|
        % kh is an array of horizontal StringState k vectors on the sphere.
        %     i.e. same latitudes, +/- longitudes.
        % ikh is the corresponding antisymmetric vectors of kh
        x, y, z, R, la, latM, longM, kv, ikv, kh, ikh, Av, Ah
    end
    
    methods
        function fs = FuzzySphere(R, N, la, nlat)
            % CONSTRUCTOR Initialize the 3 NxN matrices
            % in the basis where fs.z is diagonal.
            % la is optional and can be a boolean or of type FSLaplacian
            % if la is an FSLaplacian then simply assign it to the la field.
            % if la is true then calculate the FSLaplacian for the sphere,
            % otherwise do not assign the la field.
            % if nlat is not zero, calculate an array of k vectos for
            % drawing/animating with resolution (nlat, nlong=2*nlat).
            if nargin > 0
                j = (N - 1) / 2;
                m = j : -1 : -j;

                A = R / sqrt(j * (j + 1));

                cp = sqrt((j - m(2:end)) .* (j + m(2:end) + 1));
                Jplus = diag(cp, 1);

                cm = sqrt((j + m(1:end-1)) .* (j - m(1:end-1) + 1));
                Jminus = diag(cm, -1);

                fs.x = A * (Jplus + Jminus) / 2;
                fs.y = A * (Jplus - Jminus) / 2i;
                fs.z = A * diag(m);
                fs.R = R;
                
                if nargin > 2 
                    if isa(la, 'FSLaplacian')
                        fs.la = la;
                    elseif la
                        fs.la = FSLaplacian(fs);
                    end
                end
                
                if nargin == 4
                    if ~(nlat == 0) && ~isempty(fs.la)
                        fs.assign_karrays(nlat)
                    end
                end
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
        
        function result = Jx(self)
            N = size(self.x, 1);
            result = sqrt(N^2-1) / 2 / self.R * self.x;
        end
        
        function result = Jy(self)
            N = size(self.x, 1);
            result = sqrt(N^2-1) / 2 / self.R * self.y;
        end
        
        function result = Jz(self)
            N = size(self.x, 1);
            result = sqrt(N^2-1) / 2 / self.R * self.z;
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
        
        function obj = assign_karrays(self, nlat)
            N = size(self.x, 1);
            nlong = 2*nlat;
            lat = linspace(-90, 90, nlat);
            long = linspace(-180, 180, nlong);
            [self.longM, self.latM] = meshgrid(long, lat);
           
            self.kv = zeros(nlat, nlong, N^2);
            self.ikv = zeros(nlat, nlong, N^2);
            self.kh = zeros(nlat, nlong, N^2);
            self.ikh = zeros(nlat, nlong, N^2);
            self.Av = zeros(nlat, nlong, N^2);
            self.Ah = zeros(nlat, nlong, N^2);
            
            % possible 2 times speed up here by using symmetry
            for ii = 1:nlat
                theta = deg2rad(lat(ii));
                for jj = 1:nlong
                    phi = deg2rad(long(jj));
                    
                    ssv = StringState(theta, phi, self);
                    
                    cs1 = CoherentState([phi, theta, 1], self);
                    cs2 = CoherentState([-phi, theta, 1], self);
                    ssh = StringState(cs1, cs2, self);
                    
                    self.kv(ii, jj, :) = ssv.k0;
                    self.ikv(ii, jj, :) = ssv.ik0;
                    self.kh(ii, jj, :) = ssh.k0;
                    self.ikh(ii, jj, :) = ssh.ik0;
                    
                    Mv0 = StringState.k2M(ssv.k0, self.la);
                    iMv0 = StringState.k2M(ssv.ik0, self.la);
                    Av0 = 0.5 * (Mv0 - 1i*iMv0);
                    Mh0 = StringState.k2M(ssh.k0, self.la);
                    iMh0 = StringState.k2M(ssh.ik0, self.la);
                    Ah0 = 0.5 * (Mh0 - 1i*iMh0);
                    
                    self.Av(ii,jj,:) = Av0(:);
                    self.Ah(ii,jj,:) = Ah0(:);
                end
            end
        end
    end
end