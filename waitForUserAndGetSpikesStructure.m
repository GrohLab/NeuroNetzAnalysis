function spikes = waitForUserAndGetSpikesStructure()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
h = gcf;
spikes = [];
set(h,'CloseRequestFcn',{@getSpikeStructureCls,h})
try
    while strcmp(h.BeingDeleted,'off')
        waitforbuttonpress
    end
catch
    disp('Figure closed. Search for the spikes')
end
end

function getSpikeStructureCls(varargin)
try
    h = varargin{end};
    figdata = get(h,'UserData');
    assignin('caller','spikes',figdata.spikes);
catch
    disp('There was a problem reading the spike structure')
    disp('Please start an issue process in GitHub...')
    disp('Sorry :(')
end
delete(h)
end