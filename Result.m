classdef Result
    % Result class to store and visualize the result of the numerical
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
        
        function frame = make_frame(obj, J, tn)
            n = 20;
            a = linspace(-pi/2, pi/2, n);
            % [A1, A2] = meshgrid(a, a);
            
            N = sqrt(size(obj.y, 2) / 2);
            A = zeros(n);
            
            for i = 1:n
                for j = 1:n
                    n1 = [0, a(i), 1];
                    n2 = [0, a(j), 1];
                    
                    cs1 = CoherentState(J, n1, CoordType.spherical);
                    cs2 = CoherentState(J, n2, CoordType.spherical);
                    Phi = reshape(obj.y(tn, 1:N^2), N, N);
                    A(i, j) = abs(ctranspose(cs1.v) * Phi * cs2.v);
                end
            end
            
            ax = gca;
            pcolor(A);
            shading interp;
            colorbar;
            set(ax, 'XTickLabel', linspace(-0.4, 0.5, 10));
            set(ax, 'YTickLabel', linspace(-0.4, 0.5, 10));
            xlabel('Angle 1 (\theta_1 / \pi)', 'FontSize', 16)
            ylabel('Angle 2 (\theta_2 / \pi)', 'FontSize', 16)
            
            frame = getframe;
        end
        
        function ani = animate(obj)
            % Animate the result
            ani = obj.Property1 + inputArg;
        end
    end
end

