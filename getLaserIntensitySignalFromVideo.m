function lsrInt = getLaserIntensitySignalFromVideo( video_obj, dlc_table, varargin )


defaultMemOcc = 1/4;
checkMemOcc = @(x) isscalar(x) && isnumeric(x) && x > 0 && x < 0.5;

p = inputParser;

addRequired(p, "video_obj", @(x) isa(x, 'VideoReader' ) );
addRequired(p, "dlc_table", @(x) isa(x, 'table' ) );
addParameter(p, "mem_occ", defaultMemOcc, checkMemOcc );

parse(p, video_obj, dlc_table, varargin{:} );

video_obj = p.Results.video_obj;
dlc_table = p.Results.dlc_table;
mem_occ = p.Results.mem_occ;

nonzero_table = @(x) nnz( x(:,1) );
Nf = cellfun(@(x) nonzero_table( dlc_table.(x) ), ...
    dlc_table.Properties.VariableNames );
if std(Nf)
    Nf = min( Nf );
else
    Nf = mean( Nf );
end
wd = video_obj.Width;
ht = video_obj.Height;
fr = video_obj.FrameRate;

lsrInt = zeros( Nf, 1, 'single' );

frameByte = ht * wd * 3;
try
    [~, mem] = memory; maxMem = mem.PhysicalMemory.Available;
    bufferFrames = floor( (maxMem/frameByte) * mem_occ );
catch
    % No memory module installed. Considering 32 GB of RAM memory
    bufferFrames = floor( (32e9/frameByte) * mem_occ );
end

frames = zeros(ht, wd, 3, bufferFrames, 'uint8');
[X, Y] = meshgrid( 1:wd, 1:ht );

% Dot product between normalised likelihood and positions to get the most
% certain position in the given time period
weighted_mean = @(x) x(:,[1,2])' * ( x(:,3) ./ vecnorm( x(:,3), 1) );
frCount = 0;
while frCount < Nf && hasFrame( video_obj )

    if frCount + bufferFrames > Nf
        bufferFrames = Nf - frCount;
    elseif frCount == Nf
        bufferFrames = 1;
    end
    mu_lcoords = weighted_mean( ...
        dlc_table.laser(frCount + [1, bufferFrames], :) ) ;
    lroi = 8 > sqrt( ( X - mu_lcoords(1) ).^2 + ( Y - mu_lcoords(2) ).^2 );
    [xroi, yroi] = find( lroi );
    fprintf(1, 'Reading %d to %d out of %d frames...\n', ...
        frCount + [1,bufferFrames], Nf)
    try
        frames(:,:,:,1:bufferFrames) = read( video_obj, ...
            frCount + [1, bufferFrames]);
    catch
        fprintf(1, "Error while calculating frames in video!\n")
        try
            video_obj.CurrentTime = frCount / fr;
            cf = 1;
            while hasFrame( video_obj )
                frames(:,:,:,cf) = readFrame( video_obj );
                frCount = frCount + 1; cf = cf + 1;
            end
        catch
        end
        % Nf = frCount;
        bufferFrames = cf;
        fprintf(1, "File limit: %d\n", frCount)
        % frames(:,:,:,1:bufferFrames) = read( video_obj, ...
        %     frCount + [1, bufferFrames]);
    end

    lsrInt(frCount + (1:bufferFrames)) = ...
        mean( squeeze( frames( xroi, yroi, 1, 1:bufferFrames ) ), [1, 2] );
    frames = zeros(ht, wd, 3, bufferFrames, 'uint8');
    frCount = frCount + bufferFrames;
end

end