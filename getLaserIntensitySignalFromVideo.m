function lsrInt = getLaserIntensitySignalFromVideo( video_obj, dlc_table )

Nf = video_obj.NumFrames;
wd = video_obj.Width;
ht = video_obj.Height;

lsrInt = zeros( Nf, 1 , 'uint8' );
mu_lcoords = round( mean( dlc_table.laser ) );

mem = memory; maxMem = mem.MemAvailableAllArrays;
frameByte = ht * wd * 3;
bufferFrames = floor((maxMem/frameByte) * 0.6);

frames = zeros(ht, wd, 3, bufferFrames, 'uint8');

frCount = 0;
while hasFrame( video_obj )
    if frCount + bufferFrames > Nf
        bufferFrames = Nf - frCount;
    end
    fprintf(1, 'Reading %d/%d frames...\n', frCount + bufferFrames, Nf)
    frames(:,:,:,1:bufferFrames) = read( video_obj, frCount + [1, bufferFrames]);
    %frames(:,:,[2,3],:) = [];
    %frames = squeeze(frames);
    lsrInt(frCount + (1:bufferFrames)) = ...
        frames( mu_lcoords(1), mu_lcoords(2), 1, 1:bufferFrames );
    frCount = frCount + bufferFrames;
end

end