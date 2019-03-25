function [F, V] = animateMatrixElements(t, ss)
%ANIMATEMATRIXELEMENTS Animate
%   t is a vector of the time steps
%   ss is the StringState object that is to be time evolved.
%   F is a vector of frames
%   V is a video file (.mp4)

fig1 = figure('Position', [10, 10, 900, 600], 'Name', 'Animation of matrix elements (theta = 0, N=50)');
cmin = -.1; cmax = .1;

ax1 = subplot(2, 3, 1);
im1 = image(abs(ss.getM), 'CDataMapping', 'scaled');
title('Absolute magnitude of matrix elements')
axis square;
caxis([0, cmax])
colorbar;

ax2 = subplot(2, 3, 2);
im2 = image(real(ss.getM), 'CDataMapping', 'scaled');
title('Real part of matrix elements')
axis square;
caxis([cmin, cmax])
colorbar;

ax3 = subplot(2, 3, 3);
im3 = image(imag(ss.getM), 'CDataMapping', 'scaled');
title('Imaginary part of matrix elements')
axis square;
caxis([cmin, cmax])
colorbar;

ax4 = subplot(2, 3, 4);
im4 = image(abs(ss.getiM), 'CDataMapping', 'scaled');
title('Absolute magnitude of matrix elements')
axis square;
caxis([0, cmax])
colorbar;

ax5 = subplot(2, 3, 5);
im5 = image(real(ss.getiM), 'CDataMapping', 'scaled');
title('Real part of matrix elements')
axis square;
caxis([cmin, cmax])
colorbar;

ax6 = subplot(2, 3, 6);
im6 = image(imag(ss.getiM), 'CDataMapping', 'scaled');
title('Imaginary part of matrix elements')
axis square;
caxis([cmin, cmax])
colorbar;

F(length(t)) = struct('cdata',[],'colormap',[]);
V = VideoWriter('../Videos/matrix_elements_animation.mp4', 'MPEG-4');
V.FrameRate = 7;
open(V);

%axes = [ax1, ax2, ax3, ax4, ax5, ax6];

for n = 1:length(t)
    Mt = StringState.k2M(ss.kt(t(n)), ss.fs.la);
    iMt = StringState.k2M(ss.ikt(t(n)), ss.fs.la);
    im1.CData = abs(Mt);
    im2.CData = real(Mt);
    im3.CData = imag(Mt);
    im4.CData = abs(iMt);
    im5.CData = real(iMt);
    im6.CData = imag(iMt);
    
    drawnow;
    frame = getframe(gcf);
    F(n) = frame;
    writeVideo(V, frame);
end

close(V);
end

