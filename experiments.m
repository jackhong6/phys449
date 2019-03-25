function experiments(n)
% Some scratch code to experiment with.
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
        %plot(t, overlaps1);
        plot(t, overlaps2);
        title_txt = sprintf('Overlap for $N=%d$, $v=%.2f$, $\\theta = $ %.2f $\\pi$', ... 
                N, v0, theta/pi);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | ik(0) \rangle$', 'interpreter', 'latex')
        
    case (5)
        % Check if other antipodal states give the same overlap  
        N = 50;
        load FuzzySpheres.mat fs50
        fs = fs50;
        v0 = 0;
        
        hold on
        % x axis antipodal state
        n1 = [1, 0, 0];
        n2 = [-1, 0, 0];
        cs1 = CoherentState(fs, n1, CoordType.cartesian);
        cs2 = CoherentState(fs, n2, CoordType.cartesian);
        ss = StringState(cs1, cs2, fs);
        ss.k0 = ss.getk;
        ss.dkdt0 = ss.getdkdt(v0, [1, 0, 0], CoordType.spherical);
        t = linspace(0, 1, 500);
        overlaps = zeros(length(t), 1);
        for ii = 1:length(t)
            kti = ss.kt(t(ii));
            overlaps(ii) = kti(:)' * ss.k0;
        end
        
        overlaps = overlaps ./ overlaps(1);
        subplot(2, 2, 1)
        plot(t, overlaps);
        title_txt = sprintf('Antipodal states (x axis, $N=%d$, $v=%.2f$)', ...
            N, v0);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | k(0) \rangle$', 'interpreter', 'latex')
        
        % y axis antipodal state
        n1 = [0, 1, 0];
        n2 = [0, -1, 0];
        cs1 = CoherentState(fs, n1, CoordType.cartesian);
        cs2 = CoherentState(fs, n2, CoordType.cartesian);
        ss = StringState(cs1, cs2, fs);
        ss.k0 = ss.getk;
        ss.dkdt0 = ss.getdkdt(v0, [1, 0, 0], CoordType.spherical);
        t = linspace(0, 1, 500);
        overlaps = zeros(length(t), 1);
        for ii = 1:length(t)
            kti = ss.kt(t(ii));
            overlaps(ii) = kti(:)' * ss.k0;
        end
        
        overlaps = overlaps ./ overlaps(1);
        subplot(2, 2, 2)
        plot(t, overlaps);
        title_txt = sprintf('Antipodal states (y axis, $N=%d$, $v=%.2f$)', ...
            N, v0);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | k(0) \rangle$', 'interpreter', 'latex')
        
        % z axis antipodal state (North - South poles)
        ss = StringState(pi/2, fs);
        ss.k0 = ss.getk;
        ss.dkdt0 = ss.getdkdt(v0, [1, 0, 0], CoordType.spherical);
        t = linspace(0, 1, 500);
        overlaps = zeros(length(t), 1);
        for ii = 1:length(t)
            kti = ss.kt(t(ii));
            overlaps(ii) = kti(:)' * ss.k0;
        end
        
        overlaps = overlaps ./ overlaps(1);
        subplot(2, 2, 3)
        plot(t, overlaps);
        title_txt = sprintf('Antipodal states (z axis, $N=%d$, $v=%.2f$)', ...
            N, v0);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | k(0) \rangle$', 'interpreter', 'latex')
                
        % x axis antipodal state
        n1 = [1, 1, 1];
        n2 = [-1, -1, -1];
        cs1 = CoherentState(fs, n1, CoordType.cartesian);
        cs2 = CoherentState(fs, n2, CoordType.cartesian);
        ss = StringState(cs1, cs2, fs);
        ss.k0 = ss.getk;
        ss.dkdt0 = ss.getdkdt(v0, [1, 0, 0], CoordType.spherical);
        t = linspace(0, 1, 500);
        overlaps = zeros(length(t), 1);
        for ii = 1:length(t)
            kti = ss.kt(t(ii));
            overlaps(ii) = kti(:)' * ss.k0;
        end
        
        overlaps = overlaps ./ overlaps(1);
        subplot(2, 2, 4)
        plot(t, overlaps);
        title_txt = sprintf('Antipodal states (n=+/-[1,1,1], $N=%d$, $v=%.2f$)', ...
            N, v0);
        title(title_txt,'interpreter', 'latex')
        xlabel('Time', 'interpreter', 'latex')
        ylabel('$\langle k(t) | k(0) \rangle$', 'interpreter', 'latex')
        
    case (6)
        N = 10;
        load FuzzySpheres.mat fs10
        fs = fs10;
        v0 = 0;
        
        n_theta = 100;
        n_t = 100;
        thetas = linspace(0, pi/2, n_theta);
        t = linspace(0, 5, n_t);
        
        w = sqrt(diag(fs.la.getFullD) + 1);
        k0 = zeros(n_theta, N^2);
        dkdt0 = zeros(n_theta, N^2);
        for ii = 1:n_theta
            ss = StringState(thetas(ii), fs);
            k0(ii, :) = ss.getk;
            dkdt0(ii, :) = ss.getdkdt(v0, [1, 0, 0], CoordType.spherical);
        end
        
        kti = zeros(n_theta, N^2);
        overlaps = zeros(n_theta, n_t);
        for ii = 1:n_t
            kti(ii, :) = k0 * cos(w(:)*t(ii)) + dkdt0 * (sin(w(:)*t(ii)) ./ w(:));
            overlaps(:, ii) = dot(kti, k0, 2);
        end
        
        ax = gca;
        pcolor(overlaps);
        colorbar;
        
    case (7)
        N = 10;
        load FuzzySpheres.mat fs10
        fs = fs10;
        v0 = 0;
        
        n_theta = 500;
        n_t = 500;
        tf = 30;
        thetas = linspace(0, pi/2, n_theta);
        t = linspace(0, tf, n_t);
        
        w = sqrt(diag(fs.la.getFullD) + 1);
        k0 = zeros(n_theta, N^2);
        dkdt0 = zeros(n_theta, N^2);
        for ii = 1:n_theta
            ss = StringState(thetas(ii), fs);
            k0(ii, :) = ss.getk;
            dkdt0(ii, :) = ss.getdkdt(0, [0, 0, 1], CoordType.cartesian);
        end
        
        ss = StringState(pi/4, fs);
        ss.k0 = ss.getk;
        ss.dkdt0 = ss.getdkdt(v0, [1, 0, 0], CoordType.cartesian);
        overlaps = zeros(n_theta, n_t);
        for ii = 1:n_t
            kti = ss.kt(t(ii));
            overlaps(:, ii) = k0 * kti(:);
        end
        
        overlaps = overlaps ./ max(max(overlaps));
        overlaps = abs(overlaps);
        %overlaps = overlaps.^2;
        ax = gca;
        pcolor(overlaps);
        shading flat
        colorbar;
        xticks(linspace(1, n_t, 5))
        xticklabels({'0', sprintf('%.1f', tf/4), sprintf('%.1f', tf/2), sprintf('%.1f', 3*tf/4), sprintf('%.1f', tf)})
        yticks([1, floor(n_theta/2), n_theta])
        yticklabels({'0', '45', '90'})
        xlabel('Time', 'FontSize', 16);
        ylabel('Opening angle (degrees)', 'FontSize', 16);
        title('|Overlap| of 45 deg string state with states of different opening angles (N=10, v=0)')
        
    case (8)
        % Trying to find where the string disappears to
        N = 10;
        load FuzzySpheres.mat fs10
        fs = fs10;
        v0 = 0;
        
        n_theta = 50;
        n_phi = 200;
        t = 10;
        thetas = linspace(0, pi/2, n_theta);
        phis = linspace(0, 2*pi, n_phi);
        
        w = sqrt(diag(fs.la.getFullD) + 1);
        k0 = zeros(n_theta, n_phi, N^2);
        for ii = 1:n_theta
            for jj = 1:n_phi
                ss = StringState(thetas(ii), phis(jj), fs);
                k0(ii, jj, :) = ss.calculate_k0;
            end
        end
        
        ss = StringState(pi/4, fs);
        kt = ss.kt(t);
        overlaps = zeros(n_theta, n_phi);
        for ii = 1:n_theta
            for jj = 1:n_phi
                overlaps(ii, jj) = dot(squeeze(k0(ii, jj, :)), kt);
            end
        end
      
        overlaps = abs(overlaps);
        overlaps = overlaps ./ max(max(overlaps));
        ax = gca;
        pcolor(overlaps);
        shading flat
        colorbar;
        xticks(linspace(1, n_phi, 5))
        xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2 \pi'})
        yticks([1, floor(n_theta/2), n_theta])
        yticklabels({'0', '\pi/4', '\pi/2'})
        xlabel('Azimuth angle', 'FontSize', 16);
        ylabel('Opening angle', 'FontSize', 16);
        title(sprintf('|Overlap| of 45 deg string state with states of different opening angles (N=10, v=0, t=%.1f)', t))
        
    case (9)
        % Animation of evolution of 45 degree StringState
        t = 0;
        N = 10;
        load FuzzySpheres.mat fs10
        fs = fs10;
        v0 = 0;
        
        n_theta = 20;
        n_phi = 4 * n_theta;
        %thetas = linspace(0, pi/2, n_theta);
        
        [lat, long] = meshgrid(linspace(-90, 90, 2*n_theta), linspace(-180, 180, n_phi));

        k01 = zeros(n_theta, n_phi, N^2);
        k02 = zeros(n_theta, n_phi, N^2);
        for ii = 1:n_theta
            theta = deg2rad(lat(1, ii));
            for jj = 1:n_phi
                phi = deg2rad(long(jj, 1));                
                ss = StringState(theta, phi, fs);
                k01(ii, jj, :) = ss.k0;
                k02(ii, jj, :) = ss.ik0;
            end
        end
        
        ss = StringState(0, fs);
        kt = ss.kt(t);
        Z = zeros(size(lat));
        
        for ii = 1:n_theta
            for jj = 1:n_phi
                overlap1 = 2*abs(squeeze(k01(ii, jj, :))' * kt);
                overlap2 = 2*abs(squeeze(k02(ii, jj, :))' * kt);
                overlap = overlap1^2 + overlap2^2;
                Z(jj, ii) = overlap;
                Z(jj, 2*n_theta - ii + 1) = overlap;
            end
        end

        ax = axesm('mollweid');
        gridm; framem;
        title('String state (\theta = \pi/2, N=10)')
        geoshow(ax, lat, long, Z, 'DisplayType', 'texturemap');
        plabel('PlabelLocation', 30);
        colorbar;
end