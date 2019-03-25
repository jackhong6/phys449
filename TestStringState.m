classdef TestStringState < matlab.unittest.TestCase
    %% Unit tests for CoherentState and StringState classes.
    
    properties
        abstol = 1e-12;  % tolerance for checking floating point equality
    end
    
    %% Test Methods Block
    methods (Test)
        function testConstructors(tc)
            fs_wo_la = FuzzySphere(1, 3, false);
            fs_w_la = FuzzySphere(1, 3, true);

            s1 = StringState(pi/4, fs_wo_la);
           
            n1 = [0, pi/4, 1];
            n2 = [0, -pi/4, 1];       
            a = CoherentState(n1, fs_wo_la, CoordType.spherical);
            b = CoherentState(n2, fs_wo_la, CoordType.spherical);
            s2 = StringState(a, b, fs_w_la);
            
            w = sqrt(diag(fs_w_la.la.getFullD) + 1);
            
            tc.verifyEqual(s1.p, s2.p)
            tc.verifyNotEmpty(s1.fs)
            tc.verifyNotEmpty(s2.fs)
            tc.verifyEqual(s1.m, 1)
            tc.verifyEqual(s2.m, 1)
            tc.verifyEmpty(s1.w)
            tc.verifyEqual(s2.w, w)
            tc.verifyEmpty(fs_wo_la.la)
            tc.verifyNotEmpty(fs_w_la.la)
        end

        function testNorthAndSouthPoleStates(tc)
            fs = FuzzySphere(sqrt(3)/2, 3);
            
            nN = [0, pi/2, 1];   % North pole
            nS = [0, -pi/2, 1];  % South pole

            a = CoherentState(nN, fs, CoordType.spherical);
            b = CoherentState(nS, fs, CoordType.spherical);
            Phi = StringState(a, b, fs);
            tc.verifyEqual(Phi.getM, 1/sqrt(2) * [0 0 1; 0 0 0; 1 0 0], 'AbsTol', tc.abstol)
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
        
        function testgetM(tc)
            fs = FuzzySphere(1, 3, false);
            p = 1:9;
            M = StringState.p2M(p);
            Phi = StringState(p, fs);
            tc.verifyEqual(Phi.getM, M, 'AbsTol', tc.abstol)
        end
        
        function testgetk0(tc)
            fs = FuzzySphere(1, 3, true);
            p = 1:9;
            Phi = StringState(p, fs);
            tc.verifyEqual(Phi.k0, FSLaplacian.p2kBasis(fs.la, p))
        end
        
        function testgetdkdt(tc)
            fs = FuzzySphere(1, 3, true);
            ss = StringState(pi/2, fs);
            dMdt = 1i * commutator(fs.z, ss.getM);
            dpdt = StringState.M2p(dMdt);
            dkdt = FSLaplacian.p2kBasis(fs.la, dpdt);
            tc.verifyEqual(ss.calculate_dkdt0(1, [0, 0, 1], CoordType.cartesian), dkdt)
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
        
        function testOverlap(tc)
            fs = FuzzySphere(1, 3, true);
            
            ss = StringState(pi/2, 0, fs);
            tc.verifyEqual(StringState.overlap(ss, ss), 1, 'AbsTol', tc.abstol)
            
            ss = StringState(pi/4, pi/3, fs);
            tc.verifyEqual(StringState.overlap(ss, ss), 1, 'AbsTol', tc.abstol)
            
            ss1 = StringState(pi/4, 0, fs);
            ss2 = StringState(pi/4, fs);
            
            tc.verifyEqual(StringState.overlap(ss1, ss2), 1, 'AbsTol', tc.abstol)
            
        end
        
%         function testAzimuthalAnglesForNorthSouthPoleStringState(tc)
%             fs = FuzzySphere(1, 3, true);
%             ss0 = StringState(pi/2, fs);
%             for azimuth = linspace(0, 2*pi, 12)
%                 ss = StringState(pi/2, azimuth, fs);
%                 tc.verifyEqual(ss.overlap0(ss0, 0), 1, 'AbsTol', 1e-3)
%             end
%         end
    end
end