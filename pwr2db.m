function dBsignal = pwr2db(freqSignal)
dBsignal = 20*log10(abs(freqSignal)/length(freqSignal));
end