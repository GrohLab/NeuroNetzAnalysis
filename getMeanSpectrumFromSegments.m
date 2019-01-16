function [sponSpect,IPCT] = getMeanSpectrumFromSegments(silIdx, sonogramStruct, fs)
tx = (0:length(silIdx) - 1) * (1/fs);
sponObj = StepWaveform(silIdx,fs);
sponSubs = sponObj.Triggers;
sponSpect = zeros(1,size(sonogramStruct.SpectrumImage,1),'single');
anglSpect = sponSpect;
for csw = 1:size(sponSubs,1)
    sponIdx = sonogramStruct.TimeAxis > tx(sponSubs(csw,1)) &...
        sonogramStruct.TimeAxis < tx(sponSubs(csw,2));
    sponSpect = sponSpect + mean(abs(sonogramStruct.SpectrumImage(:,sponIdx)),2)';
    % anglSpect = anglSpect + mean(angle(sonogramStruct.SpectrumImage(:,sponIdx)),2)';
    cmpxAngl = cmpxAngl + mean(exp(1i.*angle(sonogramStruct.SpectrumImage(:,sponIdx))),2)';
    % IPCT = IPCT + abs(mean(exp(1i*angle(sonogramStruct.SpectrumImage(:,sponIdx))),2))';
end
sponSpect = (sponSpect/csw) .* exp(1i * (anglSpect/csw));
IPCT = IPCT / csw;
end
