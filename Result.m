classdef Result
    % RESULT class to store and visualize the result of the numerical
    % simulation of the time evolution of a string state.
    
    properties
        t, y
    end
    
    methods
        function obj = Result(t, y)
            % Construct an instance of this class
            obj.t = t;
            obj.y = y;
        end
        
        function sph_frame = make_sph_frame(obj, J, tn)
            n = 20;
            [x, y, z] = sphere(n);
            
        end
        
        function frame = make_frame(obj, fs, tn)
            n = 10;
            a = linspace(-pi/2, pi/2, n);
            % [A1, A2] = meshgrid(a, a);
            
            N = sqrt(size(obj.y, 2) / 2);
            A = zeros(n);
            
            for i = 1:n
                for j = 1:n
                    n1 = [0, a(i), 1];
                    n2 = [0, a(j), 1];
                    
                    cs1 = CoherentState(fs, n1, CoordType.spherical);
                    cs2 = CoherentState(fs, n2, CoordType.spherical);
                    Phi = reshape(obj.y(tn, 1:N^2), N, N);
                    A(i, j) = ctranspose(cs1.v) * Phi * cs2.v;
                end
            end
            
            ax = gca;
            pcolor(A);
            shading interp;
            colorbar;
            caxis([-1, 1]);
            set(ax, 'XTickLabel', linspace(-0.4, 0.5, 10));
            set(ax, 'YTickLabel', linspace(-0.4, 0.5, 10));
            xlabel('Elevation Angle 1 (\theta_1 / \pi)', 'FontSize', 16);
            ylabel('Elevation Angle 2 (\theta_2 / \pi)', 'FontSize', 16);
            
            frame = getframe(gcf);
        end
        
        function F = animate(obj, J)
            % Animate the result
            nframes = length(obj.t);
            F(nframes) = struct('cdata',[],'colormap',[]);
            for k = 1:length(obj.t)
                F(k) = obj.make_frame(J, k);
            end
        end
    end
end

