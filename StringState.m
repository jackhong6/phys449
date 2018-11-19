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
            obj.M = a.v * ctranspose(b.v) + b.v * ctranspose(a.v);
        end
    end
end