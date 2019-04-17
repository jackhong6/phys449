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
        % ip and ik0 are the imaginary StringState
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
        p, ip, fs, k0, ik0, dkdt0, m, w
    end
    
    methods
        function obj = StringState(varargin)
            % CONSTRUCTOR Create an instance of this class using either two
            % CoherentStates, a p vector parametrization, or an opening angle.
            if nargin == 2
                obj.fs = varargin{2};
                if isscalar(varargin{1})
                    % Interpret first arg as an opening angle in radians
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
                    % Interpret first 2 args as CoherentStates
                    M = (varargin{1}.v(:) * varargin{2}.v(:)') ...
                        + (varargin{2}.v(:) * varargin{1}.v(:)');
                    
                    iM = 1i*(varargin{1}.v(:) * varargin{2}.v(:)') ...
                         - 1i*(varargin{2}.v(:) * varargin{1}.v(:)');
                    
                    % Normalize the string state so that Tr(M' * M) = 1
                    A = trace(M' * M);
                    M = M / sqrt(A);
                    iM = iM / sqrt(A);
                    
                    obj.p = obj.M2p(M);
                    obj.ip = obj.M2p(iM);
                else
                    % Interpret as (opening angle, azimuthal angle, etc.)
                    n1 = [varargin{2}, varargin{1}, 1];  % [azimuth, elevation, radius]
                    n2 = [varargin{2}, -varargin{1}, 1];
                    cs1 = CoherentState(n1, obj.fs, CoordType.spherical);
                    cs2 = CoherentState(n2, obj.fs, CoordType.spherical);
                    obj = StringState(cs1, cs2, obj.fs);
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
            
            if ~isempty(obj.fs.la) && ~isempty(obj.ip)
                obj.ik0 = obj.calculate_ik0;
            end
        end
        
        function M = getM(self)
            M = self.p2M(self.p);
        end
        
        function Mt = getMt(self, t)
            Mt = self.k2M(self.kt(t), self.fs.la);
        end
        
        function iM = getiM(self)
            iM = self.p2M(self.ip);
        end
        
        function iMt = getiMt(self, t)
            iMt = self.k2M(self.ikt(t), self.fs.la);
        end
        
        function w = getw(self)
            w = sqrt(diag(self.fs.la.getFullD) + self.m^2);
        end
        
        function k0 = calculate_k0(self)
            k0 = FSLaplacian.p2kBasis(self.fs.la, self.p);
        end
        
        function ik0 = calculate_ik0(self)
            ik0 = FSLaplacian.p2kBasis(self.fs.la, self.ip);
        end
        
        function dkdt0 = calculate_dkdt0(self, v, n, coordType)
            % Return a vector of velocities for the k vector given
            % initial speed v and axis of rotation n.
            % n is assumed to be a unit vector pointing in the direction
            % of the axis of rotation. The direction of rotation is given
            % by the right hand rule.
            dMdt0 = 1i * v * commutator(self.fs.dot(n, coordType), self.getM);
            dpdt0 = StringState.M2p(dMdt0);
            dkdt0 = FSLaplacian.p2kBasis(self.fs.la, dpdt0);
        end
        
        function k_t = kt(self, t)
            % Return the time evolution k0 at time t. t must be a scalar.
            % Must assign properties w, k0, and dkdt0 before using.
            k_t = self.k0 .* cos(self.w*t) + (self.dkdt0 ./ self.w) .* sin(self.w*t);
        end
        
        function ik_t = ikt(self, t)
            % Return the time evolution ik0 at time t. t must be a scalar.
            % Must assign properties w, k0, and dkdt0 before using.
            ik_t = self.ik0 .* cos(self.w*t) + (self.dkdt0 ./ self.w) .* sin(self.w*t);
        end
        
        function frame = draw(self, t)
            % DRAW Plot a string state using a map projection (mollweid).
            % ss StringState to be drawn
            % t is the time to draw (scalar)
            Z = zeros(size(self.fs.latM));

            for ii = 1:size(Z, 1)
                for jj = 1:size(Z, 2)
                    kt = self.kt(t);
                    overlap_v  = 2*abs(dot(squeeze(self.fs.kv(ii, jj, :)), kt));
                    overlap_iv = 2*abs(dot(squeeze(self.fs.ikv(ii, jj, :)), kt));
                    overlap_h  = 2*abs(dot(squeeze(self.fs.kh(ii, jj, :)), kt));
                    overlap_ih = 2*abs(dot(squeeze(self.fs.ikh(ii, jj, :)), kt));
                    overlap = sqrt(overlap_v^2 + overlap_iv^2 + overlap_h^2 + overlap_ih^2);
                    Z(ii, jj) = overlap;
                end
            end

            axesm('mollweid');
            %title('String state')
            geoshow(self.fs.latM, self.fs.longM, Z, 'DisplayType', 'texturemap');
            %plabel('PlabelLocation', 30);
            %colorbar;
            caxis([0, 1])
            frame = getframe(gcf);
        end
            
        function frame = drawA(self, t)
            % DRAWA Draw the non-hermitian string state A = 0.5 * (M - iM) 
            %       using a map projection (mollweid).
            % ss StringState to be drawn. 
            % t is the time to draw (scalar)
            Z = zeros(size(self.fs.latM));
                
            for ii = 1:size(Z, 1)
                for jj = 1:size(Z, 2)
                    kt = self.kt(t);
                    ikt = self.ikt(t);
                    Mt = StringState.k2M(kt, self.fs.la);
                    iMt = StringState.k2M(ikt, self.fs.la);
                    A = 0.5 * (Mt - 1i*iMt);

                    overlap_v  = abs(A(:)' * squeeze(self.fs.Av(ii, jj, :)));
                    %overlap_h = abs(A(:)' * squeeze(self.fs.Ah(ii, jj, :)));
                    overlap = sqrt(overlap_v^2);% + overlap_h^2);
                    Z(ii, jj) = overlap;
                end
            end

            axesm('mollweid');
            %title('String state')
            geoshow(self.fs.latM, self.fs.longM, Z, 'DisplayType', 'texturemap');
            %plabel('PlabelLocation', 30);
            %colorbar;
            frame = getframe(gcf);
        end
        
        function [F, V] = animate(self, t)
            % ANIMATE Create video of the string state at the given times.
            %         t is an array of time steps
            %         F is an array of frames
            %         V is the video in MP4 format
            figure();
            axesm('mollweid');
            F(length(t)) = struct('cdata',[],'colormap',[]);
            colorbar;
            
            V = VideoWriter('../Videos/animateStringState.mp4', 'MPEG-4');
            V.FrameRate = 5;
            
            open(V);
            for n = 1:length(t)
                frame = self.draw(t(n));
                F(n) = frame;
                writeVideo(V, frame);
            end
            close(V);
            
        end
        
        function [F, V] = animateA(self, t)
            % ANIMATE Create video of the string state at the given times.
            %         t is an array of time steps
            %         F is an array of frames
            %         V is the video in MP4 format
            figure();
            axesm('mollweid');
            F(length(t)) = struct('cdata',[],'colormap',[]);
            colorbar;
            
            V = VideoWriter('../Videos/animateStringStateA.mp4', 'MPEG-4');
            V.FrameRate = 5;
            
            open(V);
            for n = 1:length(t)
                frame = self.drawA(t(n));
                F(n) = frame;
                writeVideo(V, frame);
            end
            close(V);
            
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
        
        function M = k2M(k, la)
            M = StringState.p2M(FSLaplacian.k2pBasis(la, k));
        end
        
        function result = overlap(varargin)
            % Return the normalized overlap between states a and b.
            if nargin == 3
                overlap1 = 2*abs(varargin{1}(:)' * varargin{2}(:));
                overlap2 = 2*abs(varargin{1}(:)' * varargin{3}(:));
                result = (overlap1^2 + overlap2^2);
            elseif isa(varargin{1}, 'StringState')
                a = varargin{1};
                b = varargin{2};
                result = StringState.overlap(a.k0, b.k0, b.ik0);
            else
                
            end
        end
    end
end