function plv=sPLV(X,fs,fmin,fmax)
 
%% Inputs/Outputs 
% X:  numChannels x NSamples --> signals
% sf: Sampling frequency
% fmin,fmax: frequency band range
% plv: N x numChannels x numChannels --> N is the number of windows
 
%% Filtering
FrequencyBand =[fmin fmax];
nchannels = size(X,1);
filterorder=6/FrequencyBand(1);
 
% filtering with FIR
b1 = fir1(floor(filterorder*fs),[FrequencyBand(1) FrequencyBand(2)]/(fs/2));
for k = 1:nchannels
    FiltSig(k,:) = filtfilt(b1,1,double(X(k,:)));
end
 
%% Extract the instantaneous phase using Hilbert Transform
numSamples= size(FiltSig,2);
f_interest=fmin+(fmax-fmin)/2;
%window = 6/(f_interest); % Use 6 cycles as recomended by Lachaux et al.2000
winSamples = numSamples; %floor(window * fs); 
N=floor(numSamples/winSamples);
plv = zeros(N,nchannels, nchannels);
 
for channelCount = 1:nchannels
    FilteredSignal(channelCount, :) = angle(hilbert((FiltSig(channelCount, :))));
end
 
%% Compute PLV
 
No=winSamples/2+1;
 
for count=1:N
    for channelCount = 1:nchannels-1
        channelData = squeeze(FilteredSignal(channelCount, No-winSamples/2:No+winSamples/2-1));
        for compareChannelCount = channelCount+1:nchannels
            compareChannelData = squeeze(FilteredSignal(compareChannelCount, No-winSamples/2:(No+winSamples/2-1)));
                diff=channelData(:, :) - compareChannelData(:, :);
                diff=diff';
               plv(count,channelCount,compareChannelCount) =abs(sum(exp(1i*diff)))/length(diff);  
        end
    end
    
    No=No+winSamples;
end
