function [iOk] = createEmptyRecordingFiles( beh_path )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnOpts = {'UniformOutput', false};
expandPath = @(x) fullfile( x.folder, x.name );

if exist( beh_path, "dir")
    tf_files = dir( fullfile( beh_path, "TriggerSignals*.bin" ) );
    if ~isempty( tf_files )
        tf_dates = getDates( tf_files, "TriggerSignals", ".bin");
        rf_names = "Recording" + string( tf_dates ) + ".bin";
        fID = arrayfun(@(x) fopen( fullfile( beh_path, x), "w" ), rf_names );
        
    end
else
end


end