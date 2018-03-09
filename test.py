#!/usr/bin/python
# -*- coding: utf-8 -*-

import scipy.io as sio
import scipy.linalg as lag
import scipy.signal as sig
import matplotlib.pyplot as plt
import numpy as np

# Memory and relevance improved correlation function:
def correlate_time_window(a, v, time_window):
    tw = time_window
    corr_signal = np.zeros(2*tw + 1)
    for k in range(-tw,tw+1):
        if k < 0:
            temp = np.dot(a[:k],v[-k:])
        else:
            if k == 0:
                temp = np.dot(a,v)
            else:
                temp = np.dot(a[k:],v[:-k])
        corr_signal[k+tw]=temp        
    return corr_signal

def time_freq_analysis(x, nperseg, noverlap):
    winL = nperseg - noverlap
    Nwin = int(np.floor(len(x)/winL))
    tfimg_m = np.zeros((nperseg, Nwin))
    tfimg_p = np.array(tfimg_m)
    for cW in range(Nwin):
        seg = x[winL*cW:winL*(cW+1)]
        seg_sp = np.fft.fftshift(np.fft.fft(seg))
        tfimg_m[cW,:] = np.absolute(seg_sp)
        tfimg_p[cW,:] = np.unwrap(np.angle(seg_sp))
    return (tfimg_p, tfimg_m)

_plot = False
# Importing the time series from the .mat file
traces = sio.loadmat('..\ForEmilio.mat')
LFP = np.array(traces['LFP'])
sp1 = np.array(traces['sp1'])
sp2 = np.array(traces['sp2'])
N = len(LFP)
LFP = LFP.reshape(N)
sp1 = sp1.reshape(sp1.size)
sp2 = sp2.reshape(sp2.size)
fs = 20e3;                        # Sampling frequency
tx = np.linspace(0.0,(N - 1)/fs,N)  # Time axis

# Decompressing the spikes: Create time series from the spike times
# ('sp1' and 'sp2)
st1 = np.zeros(N)
st2 = np.zeros(N)
st1[sp1-1] = 1
st2[sp2-1] = 1
POm_t = sp1/fs
L5_t = sp2/fs

# Correlation of 10000 milliseconds forth and back for two direction
# causality from the both spike trains.

# The result is around 8 time frames, which means 0.4 ms delay for the
# synapse. This result gives us only 13 synapses effectively transmitted
# (mx_corrcoef = 13).
# URL:http://www.oxfordreference.com/view/10.1093/oi/authority.20110803100547277

time_window = 1
lag_samples = int(time_window*fs)
lgx = np.array(range(-lag_samples,lag_samples+1))
corr_signal = correlate_time_window(st1,st2,lag_samples)
mx_corrcoef = corr_signal.max()
mx_lag = lgx[corr_signal==mx_corrcoef]
eff_syn = mx_lag/fs #st1-->st2
if _plot:
    lgnd = "(" + str(mx_lag)+", " + str(mx_corrcoef)+")"
    plt.plot(lgx,corr_signal);plt.text(mx_lag,mx_corrcoef,
                                       lgnd);
    plt.xlabel('Lag [k]');plt.ylabel('Correlation coefficient')
    plt.title(r"Cross-correlation $(POm\star L5B)$")
    plt.grid(True);plt.show()

st1_n = st1*(1/np.sqrt(st1.sum()))
st2_n = st2*(1/np.sqrt(st2.sum()))
LFP_n = LFP*(1/lag.norm(LFP))

corr_field1 = correlate_time_window(-st1_n,LFP_n,lag_samples)
mx_cc_field1 = corr_field1.max()
mx_l_field1 = lgx[corr_field1==mx_cc_field1]
lgnd1 = "(" + str(mx_l_field1)+", " + str(mx_cc_field1)+")"
if _plot:
    plt.plot(lgx,corr_field1);plt.text(mx_l_field1,mx_cc_field1,
                                       lgnd1)
    plt.title(r"Cross-correlation $(POm\star LFP)$")
    plt.ylabel("Normalized correlation coefficients")
    plt.xlabel("Lag [k]");plt.grid(True);plt.show()

corr_field2 = correlate_time_window(-st2_n,LFP_n,lag_samples)
mx_cc_field2 = corr_field2.max()
mx_l_field2 = lgx[corr_field2 == mx_cc_field2]
if _plot:
    lgnd2 = "(" + str(mx_l_field2)+", " + str(mx_cc_field2) + ")"
    plt.plot(lgx,corr_field2);plt.text(mx_l_field2,mx_cc_field2,
                                       lgnd2);
    plt.title(r"Cross-correlation $(L5B\star LFP)$")
    plt.ylabel("Normalized correlation coefficients")
    plt.xlabel("Lag [k]");plt.grid(True);plt.show()

# Frequency analysis of the signals. Searching for oscillatory
# behaviours

##fx = np.fft.fftshift(np.fft.fftfreq(N, 1/fs))
##LFPs = np.fft.fftshift(np.fft.fft(LFP))
# Coherence between the different signals. Does it makes sense?

# Zero padding for faster performance of the FFT.
zp = int(2**np.math.ceil(np.log2(st1_n.size))) - st1_n.size
# Frequency axis and zero padded signals for faster performance.
fx = np.fft.fftshift(np.fft.fftfreq(N + zp,1/fs))
st1_zp = np.pad(st1_n,(0,zp),'constant')
st2_zp = np.pad(st2_n,(0,zp),'constant')
LFP_zp = np.pad(LFP_n,(0,zp),'reflect')#<- To avoid introducing high
# frequency noise

# Signal spectrum
st1_sp = np.fft.fftshift(np.fft.fft(st1_zp))
st2_sp = np.fft.fftshift(np.fft.fft(st2_zp))
LFP_sp = np.fft.fftshift(np.fft.fft(LFP_zp))

# Spectrum's magnitude (results contain both the frequency axis and the
# power spectrum)
st1_ps = sig.periodogram(st1_zp,fs,window='hann')
st2_ps = sig.periodogram(st2_zp,fs,window='hann')
LFP_ps = sig.periodogram(LFP_zp,fs,window='hann')

st1_m = 20 * np.log10(np.absolute(st1_sp))
st2_m = 20 * np.log10(np.absolute(st2_sp))
LFP_m = 20 * np.log10(np.absolute(LFP_sp))
# Spectrum's phase
st1_p = np.unwrap(np.angle(st1_sp))
st2_p = np.unwrap(np.angle(st2_sp))
LFP_p = np.unwrap(np.angle(LFP_sp))





fx_s, Cpl = sig.coherence(st1_zp,st2_zp,fs,window='hann',
                        nperseg=1024,noverlap=256)
Cplf = sig.coherence(st1_zp,LFP_zp,fs,window='hann',
                        nperseg=1024,noverlap=256)
Cllf = sig.coherence(st2_zp,LFP_zp,fs,window='hann',
                        nperseg=1024,noverlap=256)



# Performing the np.welch method to the signals, I can see some patterns
# that might be quite useful. Consider this. Explore the phase as well.

