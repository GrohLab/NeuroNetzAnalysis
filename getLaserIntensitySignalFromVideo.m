function lsrInt = getLaserIntensitySignalFromVideo( video_obj, dlc_table )

Nf = video_obj.NumFrames;
wd = video_obj.Width;
ht = video_obj.Height;

lsrInt = zeros( Nf, 1, 'single' );

frameByte = ht * wd * 3;

try
    [~, mem] = memory; maxMem = mem.PhysicalMemory.Available;
    bufferFrames = floor( (maxMem/frameByte) * 0.8 );
catch
    % No memory module installed. Considering 32 GB of RAM memory
    bufferFrames = floor( (32e9/frameByte) * 0.5 );
end

frames = zeros(ht, wd, 3, bufferFrames, 'uint8');
[X, Y] = meshgrid( 1:wd, 1:ht );

frCount = 0;
while frCount < Nf
    
    if frCount + bufferFrames > Nf
        bufferFrames = Nf - frCount;
    elseif frCount == Nf
        bufferFrames = 1;
    end
    mu_lcoords = mean( dlc_table.laser(frCount + [1, bufferFrames], [1,2]) );
    lroi = 8 > sqrt( ( X - mu_lcoords(1) ).^2 + ( Y - mu_lcoords(2) ).^2 );
    [xroi, yroi] = find( lroi );
    fprintf(1, 'Reading %d-%d out of %d frames...\n', ...
        frCount + [1,bufferFrames], Nf)
    frames(:,:,:,1:bufferFrames) = read( video_obj, ...
        frCount + [1, bufferFrames]);
    %frames(:,:,[2,3],:) = [];
    %frames = squeeze(frames);
    lsrInt(frCount + (1:bufferFrames)) = ...
        mean( squeeze( frames( xroi, yroi, 1, 1:bufferFrames ) ), [1, 2] );
    frCount = frCount + bufferFrames;
end

end