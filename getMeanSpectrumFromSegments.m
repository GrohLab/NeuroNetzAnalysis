function [sponSpect,ITPC] = getMeanSpectrumFromSegments(silIdx, sonogramStruct, fs)
sfs = 1/mean(diff(sonogramStruct.TimeAxis));
% tx = (0:length(silIdx) - 1) * (1/fs);
sponObj = StepWaveform(silIdx,fs);
sponSubs = sponObj.Triggers;
sSubs = round(sponSubs * (sfs/fs));
sSubs(sSubs == 0) = 1;
Nt = max(diff(sSubs,1,2))+1;
sponSpect = zeros(... Rows --> Frequency, columns --> time
    size(sonogramStruct.SpectrumImage,1),Nt,'single');
% vectSpect = sponSpect;
% partialITPC = anglSpect;
cmpxAngl = sponSpect;
for csw = 1:size(sponSubs,1)
%     sponIdx = sonogramStruct.TimeAxis >= tx(sponSubs(csw,1)) &...
%         sonogramStruct.TimeAxis <= tx(sponSubs(csw,2));
%     if sum(sponIdx) ~= Nt
     sRng = sSubs(csw,1):sSubs(csw,2);
     if diff(sSubs(csw,:))+1 ~= Nt
        if csw < size(sponSubs,1) % If it is the beguinning of the signal
            rRng = Nt-diff(sSubs(csw,:)):Nt;
            sponSpect(:,rRng) = sponSpect(:,rRng) +...
                abs(sonogramStruct.SpectrumImage(:,sRng));
            cmpxAngl(:,rRng) = cmpxAngl(:,rRng) +...
                exp(1i.*angle(sonogramStruct.SpectrumImage(:,sRng)));
        else % If it is the end (most probably)
            rRng = 1:diff(sSubs(csw,:))+1;
            sponSpect(:,rRng) = sponSpect(:,rRng) +...
                abs(sonogramStruct.SpectrumImage(:,sRng));
            cmpxAngl(:,rRng) = cmpxAngl(:,rRng) +...
                exp(1i.*angle(sonogramStruct.SpectrumImage(:,sRng)));
        end
    else
        sponSpect = sponSpect +...
            abs(sonogramStruct.SpectrumImage(:,sRng));
        cmpxAngl = cmpxAngl +...
            exp(1i.*angle(sonogramStruct.SpectrumImage(:,sRng)));
    end
end
sponSpect = (sponSpect./csw) .* exp(1i * (angle(cmpxAngl./csw)));
ITPC = abs(cmpxAngl./csw);
end
