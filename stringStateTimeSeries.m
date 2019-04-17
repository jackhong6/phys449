ii = 1;
for t = linspace(0, 6, 9)
    subaxis(9, 1, ii, 'm', 0, 'p', 0, 's', 0, 'holdaxis', 1);
    ss0.draw(t)
    ii = ii + 1
%     subaxis(9, 2, ii+1, 'm', 0, 'p', 0, 's', 0, 'holdaxis', 1);
%     ss45h.draw(t)
%     ii = ii + 2;
end