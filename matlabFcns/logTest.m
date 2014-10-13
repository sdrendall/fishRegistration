
[y, x] = size(nseg);
[xx, yy] = meshgrid(1:x, 1:y);

slantSeg = nseg + xx*(.5./x);

figure, imshow(slantSeg, [])

figure
%for sig = 1:.5:10
sig = 5;
    f = fspecial('log', 150, sig);
    fseg = imfilter(nseg, f, 'symmetric');
    subplot(1,2,1)
    imagesc(fseg)
    colorbar
    title(['sigma = ', num2str(sig)])

    fsegSlant = imfilter(slantSeg, f, 'symmetric');
    subplot(1,2,2)
    imagesc(fsegSlant)
    colorbar

    %figure
    %imagesc(fseg - fsegSlant)
    %colorbar
    pause
%end