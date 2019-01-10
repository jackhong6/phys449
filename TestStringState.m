classdef TestStringState < matlab.unittest.TestCase
    %% Unit tests for CoherentState and StringState classes.
    
    properties
        abs_tol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function testNorthAndSouthPoleStates(tc)
            fs = FuzzySphere(sqrt(3)/2, 3);
            
            nN = [0, pi/2, 1];   % North pole
            nS = [0, -pi/2, 1];  % South pole

            a = CoherentState(fs, nN, CoordType.spherical);
            b = CoherentState(fs, nS, CoordType.spherical);
            Phi = StringState(a, b);
            tc.verifyEqual(a.v, [1; 0; 0], 'AbsTol', tc.abs_tol)
            tc.verifyEqual(b.v, [0; 0; 1], 'AbsTol', tc.abs_tol)
            tc.verifyEqual(Phi.getM, 1/sqrt(2) * [0 0 1; 0 0 0; 1 0 0], 'AbsTol', tc.abs_tol)
        end
        
        function testp2M(tc)
            p = 1:4;
            M = StringState.p2M(p);
            tc.verifyEqual(M, [sqrt(2) 3+4i; 3-4i 2*sqrt(2)])
            
            
            p = 1:9;
            M = StringState.p2M(p);
            tc.verifyEqual(M, [sqrt(2) 4+5i 8+9i; 
                               4-5i 2*sqrt(2) 6+7i; 
                               8-9i 6-7i 3*sqrt(2)])
        end
        
        function testp2MRandom(tc)
            p = unifrnd(-ones(16, 1), ones(16, 1)); % 4x1 random vector
            M = StringState.p2M(p);
            tc.verifyEqual(M, [sqrt(2)*p(1)  , p(5)+1i*p(6)  , p(11)+1i*p(12), p(15)+1i*p(16);
                               p(5)-1i*p(6)  , sqrt(2)*p(2)  , p(7)+1i*p(8)  , p(13)+1i*p(14);
                               p(11)-1i*p(12), p(7)-1i*p(8)  , sqrt(2)*p(3)  , p(9)+1i*p(10);
                               p(15)-1i*p(16), p(13)-1i*p(14), p(9)-1i*p(10) , sqrt(2)*p(4)])
        end
        
        function testM2p(tc)
            M = [sqrt(2) 3+4i; 3-4i 2*sqrt(2)];
            p = StringState.M2p(M);
            tc.verifyEqual(p', 1:4);
            
            M = [sqrt(2) 4+5i 8+9i; 
                 4-5i 2*sqrt(2) 6+7i; 
                 8-9i 6-7i 3*sqrt(2)];
            p = StringState.M2p(M);
            tc.verifyEqual(p', 1:9);
        end
        
        function testM2pRandom(tc)
            M = rand(4) + 1i * rand(4);
            M = (M + M') / 2;
            p = StringState.M2p(M);
            tc.verifyEqual(p, [diag(M)/sqrt(2); 
                real(M(1,2));
                imag(M(1,2));
                real(M(2,3));
                imag(M(2,3));
                real(M(3,4));
                imag(M(3,4));
                real(M(1,3));
                imag(M(1,3));
                real(M(2,4));
                imag(M(2,4));
                real(M(1,4));
                imag(M(1,4))
                ])
            tc.verifyEqual(p, [diag(M)/sqrt(2);
                real(M(2,1));
                -imag(M(2,1));
                real(M(3,2));
                -imag(M(3,2));
                real(M(4,3));
                -imag(M(4,3));
                real(M(3,1));
                -imag(M(3,1));
                real(M(4,2));
                -imag(M(4,2));
                real(M(4,1));
                -imag(M(4,1))
                ])
        end
    end
end