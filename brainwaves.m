function [varargout] = brainwaves(EEGsignals,sampling_frequency,varargin)
%% HELP BRAINWAVES
% Returns the requested brain waves (alpha, beta, gamma, delta and theta)
% from the given EEG signal(s). With predefined cut frequencies 
% 
% SYNTHAXIS AND DESCRIPTION
% [alpha, beta, gamma, delta, theta] = brainwaves(EEGsignals,
% sampling_frequency)
% [alpha, beta, theta] = brainwaves(EEGsignals, sampling_frequency,
% {'alpha'},'beta','theta')
% [alpha, theta] = brainwaves(EEGsignals, sampling_frequency,
% {'alpha',6,10},'theta')
% 
% This function returns the requested brainwaves into the different output 
% arguments. The used filters are two 3 order butterworth with cutoff
% frequency predifined:
% Alpha: 7  - 13 Hz
% Beta:  13 - 30 Hz 
% Gamma: 30 - 70 Hz 
% Delta: 1  - 4  Hz 
% Theta: 4  - 7  Hz
% 
% [alpha, beta, gamma, delta, theta] = brainwaves(EEGsignals,
% sampling_frequency) returns in each output argument
% (alpha, beta, ... , theta) a n-signal array for the input n-EEGsignals
% with the default sample_frequency. The output of the function is always 
% in the order of the greek alphabet (alpha, beta, gamma, delta, theta). 
% Please notice that EEGsignals must be a signal array. This means that the
% rows must be the signals.
% 
% [alpha, beta, theta] = brainwaves(EEGsignals, sampling_frequency,
% {'alpha'},'beta','theta') returns only the desired brainwaves with
% default cutoff frequency. The functions skips the waves that are not
% selected.
% 
% [beta, theta] = brainwaves(EEGsignals, sampling_frequency,
% {'beta',Wn1,Wn2},'theta') returns the selected brainwaves with customized
% cutoff frequencies defined by Wn1 and Wn2, where Wn1 is the lower cutoff
% frequency of the passband and Wn2 the upper limit. The frequencies of the
% other brainwaves are not affected.
% 
% Please notice that the Nyquist sampling theorem should be fulfilled when
% inputing the high-cut-off frequencies.

%% Input / Output arguments validation and information extraction
% Checking the number of input and output arguments
% narginchk(2,8)
% nargoutchk(1,5)
% Information extraction from the input arguments (waves and frequency
% limits
wavename = {'alpha', 'beta', 'gamma', 'delta', 'theta'};
% The upper row contains the upper frequency limits
%          Al Be Ga De Th 
cutfreq = [13 30 70 4  7;...
           7  13 30 0.5 4];
EEGchannels = size(varargin,2);
if EEGchannels > 5
    error('Unexpected input waves');
elseif EEGchannels == 0
    waveflag = ones(1,5);
else
    waveflag = zeros(1,5); % (alpha,beta,gamma,delta,theta)
    if EEGchannels > nargout
        disp('Warning: Not enough output arguments')
    elseif EEGchannels < nargout
        disp('Warning: Some output arguments might not be used');
    end
    for k=1:EEGchannels
        %     Wavename extraction
        datatype = class(varargin{k});
        switch datatype
            case 'cell'
                if size(varargin{k},2) == 1
                    deforfix = 'def';
                elseif size(varargin{k},2) == 3
                    deforfix = 'fix';
                else
                    error('Wave cut frequency error')
                end
                wavestrg = cell2mat(varargin{k}(1));
            case 'char'
                deforfix = 'def';
                wavestrg = varargin{k};
            otherwise
                error('Unexpected input: %d',varargin{k});
        end
        %     Definig lower and upper limits for the bandpass filter
        switch wavestrg
            case 'alpha'
                waveflag(1) = waveflag(1)+1;
                if waveflag(1) > 1
                    error('Alpha wave redifined');
                end
                if strcmp(deforfix,'def')
                    cutfreq(2,1) = 7;
                    cutfreq(1,1) = 13;
                else
                    cutfreq(2,1) = cell2mat(varargin{k}(2));
                    cutfreq(1,1) = cell2mat(varargin{k}(3));
                    if cutfreq(2,1) >= cutfreq(1,1)
                        error('Alpha cut frquencys are inverted or equal')
                    elseif cutfreq(2,1) <= 0 || cutfreq(1,1) <= 0
                        error('Negative cut frequency')
                    end
                end
            case 'beta'
                waveflag(2) = waveflag(2)+1;
                if waveflag(2) > 1
                    error('Beta wave redifined');end
                if strcmp(deforfix,'def')
                    cutfreq(2,2) = 13;
                    cutfreq(1,2) = 30;
                else
                    cutfreq(2,2) = cell2mat(varargin{k}(2));
                    cutfreq(1,2) = cell2mat(varargin{k}(3));
                    if cutfreq(2,2) >= cutfreq(1,2)
                        error('Beta cut frquencys are inverted or equal')
                    elseif cutfreq(2,2) <= 0 || cutfreq(1,2) <= 0
                        error('Negative cut frequency')
                    end
                end
            case 'gamma'
                waveflag(3) = waveflag(3)+1;
                if waveflag(3) > 1
                    error('Gamma wave redifined');end
                if strcmp(deforfix,'def')
                    cutfreq(2,3) = 30;
                    cutfreq(1,3) = 100;
                else
                    cutfreq(2,3) = cell2mat(varargin{k}(2));
                    cutfreq(1,3) = cell2mat(varargin{k}(3));
                    if cutfreq(2,3) >= cutfreq(1,3)
                        error('Gamma cut frquencys are inverted or equal')
                    elseif cutfreq(2,3) <= 0 || cutfreq(1,3) <= 0
                        error('Negative cut frequency')
                    end
                end
            case 'delta'
                waveflag(4) = waveflag(4)+1;
                if waveflag(4) > 1
                    error('Delta wave redifined');end
                if strcmp(deforfix,'def')
                    cutfreq(2,4) = 0.01;
                    cutfreq(1,4) = 3.5;
                else
                    cutfreq(2,4) = cell2mat(varargin{k}(2));
                    cutfreq(1,4) = cell2mat(varargin{k}(3));
                    if cutfreq(2,4) >= cutfreq(1,4)
                        error('Delta cut frquencys are inverted or equal')
                    elseif cutfreq(2,4) <= 0 || cutfreq(1,4) <= 0
                        error('Negative cut frequency')
                    end
                end
            case 'theta'
                waveflag(5) = waveflag(5)+1;
                if waveflag(5) > 1
                    error('Theta wave redifined');end
                if strcmp(deforfix,'def')
                    cutfreq(2,5) = 4;
                    cutfreq(1,5) = 7;
                else
                    cutfreq(2,5) = cell2mat(varargin{k}(2));
                    cutfreq(1,5) = cell2mat(varargin{k}(3));
                    if cutfreq(2,5) >= cutfreq(1,5)
                        error('Theta cut frquencys are inverted or equal')
                    elseif cutfreq(2,5) <= 0 || cutfreq(1,5) <= 0
                        error('Negative cut frequency')
                    end
                end
            otherwise
                error('Unknown brainwave: %s', wavestrg)
        end
    end
end

%% Wave filtering
% whm [waves hiper-matrix]
[m, n] = size(EEGsignals);
whm = zeros(m,n,5);

% [b50,a50] = butter(2,[(2*49)/sampling_frequency,...
%     (2*51)/sampling_frequency],'stop');
% 50 Hz notch filter.
% The cheby filter has too many ripples in the passing band. EIC 2018
% [bhpf,ahpf] = cheby2(4,60,0.4/sampling_frequency,'high');
% % 0.4 Hz (HIGH)
% [blpf,alpf] = cheby2(2,60,100/sampling_frequency,'low');
% 100 Hz (LOW)

% w [wave number -> 1 alpha, 2 beta, 3 gamma, 4 delta and 5 theta]
for w = 1:5
    if waveflag(w)% == 1
        % Filtering signal by signal of the EEG signal array. (zero
        % shifting filter)
        [b,a] = butter(2,[(2*cutfreq(2,w))/sampling_frequency,...
            (2*cutfreq(1,w))/sampling_frequency]);        
        for c=1:m
            % aux = filtfilt(b50,a50,EEGsignals(c,:));
            % EEGsignals(c,:) = filtfilt(blpf,alpf,EEGsignals(c,:));
            % EEGsignals(c,:) = filtfilt(bhpf,ahpf,EEGsignals(c,:));
            whm(c,:,w) = filtfilt(b,a,EEGsignals(c,:));
        end
        disp(['Filtering wave: ',wavename{w}])
    end
end

%% Preparing results
% Fixing hiper matrix
w = 1;
while size(whm,3) >= w
    if sum(whm(:,:,w)) == 0
        whm(:,:,w) = [];
    else
        w = w + 1;
    end
end

% Sending resuls to each output argument
varargout = cell(nargout,1);
for k=1:nargout
    varargout(k) = {whm(:,:,k)};
end
end