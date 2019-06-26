function discreteTraces = gatherDiscreteSignals(allEvents, fs, td)
%GATHERDISCRETESIGNALS returns a logic time trace of all the events in the
%given file. This is a simple approach without validating many possible
%(and logically correct) inputs to the function
%   The input for this function are a cell array containing the time stamps
%   (idx or seconds) or the logical train for the complete experimental (or
%   recorded) time, and a time factor. This later means either the user
%   provides the time stamps and the sampling frequency or the logic train
%   and the (exact) total experimental duration.
%   The output of the function is a NxM logical matrix, where N is the
%   number of given events and M is the number of time samples.

Ne = numel(allEvents);
sLen = max(cellfun(@length,allEvents));
bulls = cellfun(@islogical,allEvents);
if sum(bulls)
    if Ne ~= sum(bulls)
        [rws,cls] = cellfun(@size,allEvents);
    end
else
    if exist('fs','var') && exist('td','var') && fs > 1
        sLen = ceil(td * fs);
    else
        disp('It is not possible to determine the real length of the signals')
        sLen = ceil(max(cellfun(@max,allEvents(~bulls)))/fs) + 1;
    end
    [rws,cls] = cellfun(@size,allEvents);
    disp('No logical train given')
end
discreteTraces = false(Ne,sLen);
for ce = 1:Ne
    if bulls(ce)
        cl = length(allEvents{ce});
        if sLen == cl
            discreteTraces(ce,:) = allEvents{ce};
        else
            discreteTraces(ce,1:cl) = allEvents{ce};
        end
    else
        if  rws(ce) == 1 || cls(ce) == 1
            % This is a time stamp vector
            discreteTraces(ce,allEvents{ce}) = true;
        else
            % This is a start-stop vector
            for ced = 1:length(allEvents{ce})
                discreteTraces(ce,...
                    allEvents{ce}(ced,1):allEvents{ce}(ced,2)) = true;
            end
        end
    end
end


end

% sum(diff(allEvents{ce},1,find([rws(ce),cls(ce)]==2,1))) == 0