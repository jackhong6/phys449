classdef StringState
    % Class that represents a string state.
    %    A string state is the outer product of two coherent states.
    
    properties
        % The matrix representation in the Jz basis.
        M
    end
    
    methods
        function obj = StringState(a, b)
            %StringState Construct an instance of this class
            %   a, b are coherent states as defined in CoherentState.m
            obj.M = (a.v * b.v') + (b.v * a.v');
            
            % Normalize the string state so that Tr(A' * A) = 1
            A = trace(obj.M' * obj.M);
            obj.M = obj.M / sqrt(A); 
        end
    end
end