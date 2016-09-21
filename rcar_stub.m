% rCAR Stub TESTING

% Futz up some data!

fprintf('\n\nFor use at work!\n\n');

filename = '702-S1-LO-19.05.2016.11.38.05.set';

lowerBound  = 2;        % These are frequencies for the filtering
upperBound  = 41;
eegChannels = 3:16;     % Emotiv EPOC/EPOC+ data channels in EDF

if regexp(filename,'set$')         
    EEG = pop_loadset(filename);
elseif regexp(filename,'edf$')
    EEG = pop_biosig(filename);
else
    error('eegdemo: File type unknown');
end

EEG_only = pop_select(EEG, 'channel', eegChannels);
EEG_only = pop_eegfilt(EEG_only, lowerBound, upperBound, [], [0], 0, 0, 'fir1', 0);

for m = 3:5
    ss      = ge_getSampleBounds(EEG, m);
    data{m} = EEG_only.data(:,ss(1):ss(2));
end

blob_pre.Fs = 128;
blob_pre.data = data{3}';    % NB: data{1} and data{2} are empty--not used!
blob_pre = ebReRefAverage(blob_pre);
blob_pre = ebReRefEPOCParietal(blob_pre);

% rCAR down here!

[dOut, carRefEst, nnCarRefEst] = robustCARreference(blob_pre.data', 1/blob_pre.Fs);