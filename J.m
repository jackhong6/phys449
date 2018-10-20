classdef JMatrices
    % Representation of NxN angular momentum operators in the Jz basis.
    %   Detailed explanation goes here
    
    properties
        % NxN matrices
        x, y, z
    end
    
    methods
        function obj = J(N)
            % Constructor. Initialize the 3 NxN J.x, J.y, and J.z matrices
            % in the J.z basis
            j = (N-1) / 2;
            m = j : -1 : -j;
            
            obj.z = diag(m);
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end
