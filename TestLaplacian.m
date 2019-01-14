classdef TestLaplacian < matlab.unittest.TestCase
    %% Unit tests for Laplacian class.
    
    properties
        abs_tol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function test2x2Case(tc)
            fs = FuzzySphere(3, 2);
            la = Laplacian(fs);
            tc.verifyEqual(full(la.K), (1/3^2) * [1 -1 0 0; -1 1 0 0; 0 0 2 0; 0 0 0 2], 'AbsTol', tc.abs_tol)
        end
        
        function test3x3Case(tc)
            fs = FuzzySphere(1, 3);
            la = Laplacian(fs);
            K1 = [2 -2 0; -2 4 -2; 0 -2 2];
            K2 = [4 0 -2 0; 0 4 0 -2; -2 0 4 0; 0 -2 0 4];
            K3 = 6 * eye(2);
            tc.verifyEqual(la.K, sparse(blkdiag(K1, K2, K3)), 'AbsTol', tc.abs_tol);
        end
    end
end