function experiments(n)

switch(n)
    case(1)
        load FuzzySpheres.mat fs50;
        fs = fs50;
        N = 50;

        % opening angle between the two coherent states
        thetas = [pi/2, pi/3, pi/4, pi/6];

        t = linspace(0, 10, 500);
        v = 0;
        v0 = v*ones(N^2, 1);

        for jj = 1:length(thetas)
            n1 = [0, thetas(jj), 1];
            n2 = [0, -thetas(jj), 1];
            cs1 = CoherentState(fs, n1, CoordType.spherical);
            cs2 = CoherentState(fs, n2, CoordType.spherical);
            ss0 = StringState(cs1, cs2, fs);
            M0 = StringState.p2M(ss0.p);
            k0 = FSLaplacian.p2kBasis(fs.la, ss0.p);
            %kt = zeros(N^2, length(t));
            overlaps = zeros(length(t), 1);
            for ti = 1:length(t)
                %kt(:, ti) = ss0.kt(t(ti), m, k0, dkdt0);
                
                kti = ss0.kt(t(ti), k0, v0);
                %overlaps(ti) = kti(:)' * k0;
                
                Mti = StringState.p2M(FSLaplacian.k2pBasis(fs.la, kti));
                %overlaps(ti) = Mti(:)' * M0(:);
                
                overlaps(ti) = cs1.v(:)' * Mti * cs2.v(:);
            end

            subplot(2, 2, jj)
            plot(t, overlaps);
            title_txt = sprintf('Overlap for $N=%d$, $v=%.2f$, $\\theta = $ %.2f $\\pi$', ... 
                N, v, thetas(jj)/pi);
            title(title_txt,'interpreter', 'latex')
            xlabel('Time', 'interpreter', 'latex')
            ylabel('$\langle cs1 | M(t) | cs2 \rangle$', 'interpreter', 'latex')
        end
        
    case (2)
        N = 3;
        fs = FuzzySphere(1, N, true);
        theta = pi/4;
        m = 1;
        t = linspace(0, 10, 100);
        v0 = 1*ones(N^2, 1);
        
        n1 = [0, theta, 1];
        n2 = [0, -theta, 1];
        cs1 = CoherentState(fs, n1, CoordType.spherical);
        cs2 = CoherentState(fs, n2, CoordType.spherical);
        ss0 = StringState(cs1, cs2, fs);
        M0 = ss0.getM;
        k0 = FSLaplacian.p2kBasis(fs.la, ss0.p);
        kt = zeros(N^2, length(t));
        overlap1 = zeros(1, length(t));
        overlap2 = zeros(1, length(t));
        overlap3 = zeros(1, length(t));
        overlap4 = zeros(1, length(t));
        overlap5 = zeros(1, length(t));

        for ti = 1:length(t)
            kti = ss0.kt(t(ti), k0, v0);
            kt(:, ti) = kti;
            overlap1(ti) = kti(:)' * k0(:);
            M = StringState.p2M(FSLaplacian.k2pBasis(fs.la, kti));
            overlap2(ti) = M(:)' * M0(:);
            overlap3(ti) = sum(sum(M .* M0));
            overlap4(ti) = trace(M' * M0);
            overlap5(ti) = sum(diag(M' * M0));
        end
        hold on
        %plot(t, kt' * k0);
        plot(t, overlap3);
        plot(t, real(overlap3));
        plot(t, imag(overlap3));
        legend('1', '2')
        
        % They are all equivalent ways of expressing the overlap between
        % two states, but have slightly different scalings. The fastest way
        % is to just calculate the inner product of the two k vectors.
    case(3)
        N = 50;
        load FuzzySpheres.mat fs50
        fs = fs50;
        
        v = 0:15;
        theta = pi/4;
        t = linspace(0, 10, 500);
        
        for jj = 1:length(v)
            v0 = v(jj)*ones(N^2, 1);
            n1 = [0, theta, 1];
            n2 = [0, -theta, 1];
            cs1 = CoherentState(fs, n1, CoordType.spherical);
            cs2 = CoherentState(fs, n2, CoordType.spherical);
            ss0 = StringState(cs1, cs2, fs);
            M0 = StringState.p2M(ss0.p);
            k0 = FSLaplacian.p2kBasis(fs.la, ss0.p);
            %kt = zeros(N^2, length(t));
            overlaps = zeros(length(t), 1);
            for ti = 1:length(t)
                %kt(:, ti) = ss0.kt(t(ti), m, k0, dkdt0);
                
                kti = ss0.kt(t(ti), k0, v0);
                %overlaps(ti) = kti(:)' * k0;
                
                Mti = StringState.p2M(FSLaplacian.k2pBasis(fs.la, kti));
                overlaps(ti) = Mti(:)' * M0(:);
                
                %overlaps(ti) = cs1.v(:)' * Mti * cs2.v(:);
            end

            subplot(4, 4, jj)
            plot(t, overlaps);
            title_txt = sprintf('Overlap for $N=%d$, $v=%.2f$, $\\theta = $ %.2f $\\pi$', ... 
                N, v(jj), theta/pi);
            title(title_txt,'interpreter', 'latex')
            xlabel('Time', 'interpreter', 'latex')
            ylabel('$\langle M(t) | M(0) \rangle_F$', 'interpreter', 'latex')
        end
        
    case (4)
        N = 10;
        load FuzzySpheres.mat fs10
        fs = fs10;
        theta = pi/2;
        v0 = 0;

        ss = StringState(theta, fs);
        k0 = ss.getk;

        n1 = [0, theta, 1];
        n2 = [0, -theta, 1];
        cs1 = CoherentState(fs, n1, CoordType.spherical);
        cs2 = CoherentState(fs, n2, CoordType.spherical);
        iM0 = 1i*(cs1.v(:) * cs2.v(:)') - 1i*(cs2.v(:) * cs1.v(:)');
        ip0 = StringState.M2p(iM0);
        ik0 = FSLaplacian.p2kBasis(fs.la, ip0);
        
        t = linspace(0, 10, 1000);
        dMdt0 = 1i * v0 * commutator(fs(3), ss.getM);
        dpdt0 = StringState.M2p(dMdt0);
        dkdt0 = FSLaplacian.p2kBasis(fs.la, dpdt0);
        overlaps1 = zeros(length(t), 1);
        overlaps2 = zeros(length(t), 1);
        kt = zeros(length(t), 1);
        for ii = 1:length(t)
            kti = ss.kt(t(ii), k0, dkdt0);
            overlaps1(ii) = kti(:)' * ik0(:);
            overlaps2(ii) = kti(:)' * k0(:);
        end
        
        hold on
        plot(t, overlaps1);
        plot(t, overlaps2);
        title_txt = sprintf('Overlap for $N=%d$, $v=%.2f$, $\\theta = $ %.2f $\\pi$', ... 
                N, v0, theta/pi);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | ik(0) \rangle$', 'interpreter', 'latex')
        
    case (5)
        
end
end