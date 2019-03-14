function [F,V] = animateMatrixProfile(t, ss)
%ANIMATEMATRIXPOFILE Create animation of the sum of the absolute magnitude 
%                    squared of each diagonal of the StringState matrix. 
%                    
%   t is a vector of times to evaluate the matrix
%   ss is a StringState
%   F is a vector of frames
%   V is a video in .mp4 format

fig = figure();
N = sqrt(length(ss.p));

matrix_profile1 = zeros(1, N);
matrix_profile2 = zeros(1, N);
diagonals = -N+1 : N-1;
for n = 1:length(diagonals)
    matrix_profile1(n) = sum(abs(diag(ss.getM, diagonals(n))).^2);
    matrix_profile2(n) = sum(abs(diag(ss.getiM, diagonals(n))).^2);
end

ax1 = subplot(2, 1, 1);
bar1 = bar(diagonals, matrix_profile1);
title('Matrix profile for |a><b| + |b><a|')
xlabel('Diagonal number')
ylabel('Sum of magnitude squared')
ylim([0, max(matrix_profile1) + 0.05])

ax2 = subplot(2, 1, 2);
bar2 = bar(diagonals, matrix_profile2);
title('Matrix profile for i|a><b| - i|b><a|')
xlabel('Diagonal number')
ylabel('Sum of magnitude squared')
ylim([0, max(matrix_profile1) + 0.05])

F(length(t)) = struct('cdata',[],'colormap',[]);
V = VideoWriter('../Videos/matrix_profile_animation.mp4', 'MPEG-4');
V.FrameRate = 7;
open(V);

for n = 1:length(t)
    Mt = StringState.k2M(ss.kt(t(n)), ss.fs.la);
    iMt = StringState.k2M(ss.ikt(t(n)), ss.fs.la);
    
    for m = 1:length(diagonals)
        matrix_profile1(m) = sum(abs(diag(Mt, diagonals(m))).^2);
        matrix_profile2(m) = sum(abs(diag(iMt, diagonals(m))).^2);
    end
    
    bar1.YData = matrix_profile1;
    bar2.YData = matrix_profile2;
    drawnow;
    frame = getframe(gcf);
    F(n) = frame;
    writeVideo(V, frame);
end

close(V);
end

