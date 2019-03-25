classdef OverlapStates
    %OVERLAPSTATES k vector of StringState to plot a StringState
    %   Class to store  an array of k vectors to overlap with a particular
    %   StringState for visualizing the state. Also store matrix of lat and
    %   long coordinates for each overlap k vector.
    
    properties
        k0, ik0, lat1, long1, lat2, long2
    end
    
    methods
        function obj = OverlapState(fs, lat1, long1, lat2, long2)
            %OVERLAPSTATE Constructor
            % fs FuzzySphere on which to construct the 
            % 
            obj.lat1 = lat1;
            k01 = zeros(n_theta, n_phi, N^2);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

