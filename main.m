%--------------------------------------------------Q %1.1-----------------------------------------------------------
% Loading data files
Q1_Data1 = load('Q1_Data1.mat');
Q1_Data2 = load('Q1_Data2.mat');
Q1_Data3 = load('Q1_Data3.mat');

% Extracting specification details for each dataset
sampling_rate_Q1_data_1 = Q1_Data1.EEG.srate % sampling rate for Q1_Data1
num_channels_Q1_data_1 = Q1_Data1.EEG.nbchan % number of channels for Q1_Data1
num_trials_Q1_data_1 = Q1_Data1.EEG.trials  % number of trials for Q1_Data1
trial_length_Q1_data_1 = size(Q1_Data1.EEG.data, 2) %  / sampling_rate1 ; % length of each trial in seconds for Q1_Data1

sampling_rate_Q1_data_2 = Q1_Data2.EEG.srate % sampling rate for Q1_Data2
num_channels_Q1_data_2 = (Q1_Data2.EEG.nbchan) % number of channels for Q1_Data2
num_trials_Q1_data_2 = (Q1_Data2.EEG.trials) % number of trials for Q1_Data2
trial_length_Q1_data_2 = size(Q1_Data2.EEG.times, 2) % / sampling_rate2; % length of each trial in seconds for Q1_Data2

sampling_rate_Q1_data_3 = Q1_Data3.EEG.srate % sampling rate for Q1_Data3
num_channels_Q1_data_3 = Q1_Data3.EEG.nbchan % number of channels for Q1_Data3
num_trials_Q1_data_3 =  Q1_Data3.EEG.trials % number of trials for Q1_Data3
trial_length_Q1_data_3 = size(Q1_Data3.EEG.times, 2) % / sampling_rate3; % length of each trial in seconds for Q1_Data3

%----------------------------------------------------Q 1.2--------------------------------------------------------------

% Extracting first channel for each dataset
Q1_Data1_channel1 = Q1_Data1.EEG.data(1, :, :); % first channel for Q1_Data1
Q1_Data2_channel1 = Q1_Data2.EEG.data(1, :, :); % first channel for Q1_Data2
Q1_Data3_channel1 = Q1_Data3.EEG.data(1, :, :); % first channel for Q1_Data3

% Finding average ERP signal for each dataset using averaging
ERP_Q1_Data1 = mean(Q1_Data1_channel1, 3); % average ERP signal for Q1_Data1
ERP_Q1_Data2 = mean(Q1_Data2_channel1, 3); % average ERP signal for Q1_Data2
ERP_Q1_Data3 = mean(Q1_Data3_channel1, 3); % average ERP signal hz2for Q1_Data3

% Seting time axis for plot
time_Q1_data_1 = 0:1/Q1_Data1.EEG.srate:(length(ERP_Q1_Data1)-1)/Q1_Data1.EEG.srate; % time axis for Q1_Data1
time_Q1_data_2 = 0:1/Q1_Data2.EEG.srate:(length(ERP_Q1_Data2)-1)/Q1_Data2.EEG.srate;
time_Q1_data_3 = 0:1/Q1_Data3.EEG.srate:(length(ERP_Q1_Data3)-1)/Q1_Data3.EEG.srate;

hz_Q1_data_1=(0:1500-1)*500/1500;
hz_Q1_data_2=(0:1600-1)*400/1600;
hz_Q1_data_3=(0:2400-1)*600/2400;

% Ploting average ERP signal in time domain and power spectrum of average ERP signal using subplot
figure
subplot(2,3,1)
plot(time_Q1_data_1, ERP_Q1_Data1)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q1\_Data1')
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
plot([0 0]+.5,get(gca,'ylim'),'k--')


subplot(2,3,4)
plot(hz_Q1_data_1, abs(fft(ERP_Q1_Data1)),"linew",2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectrum of Average ERP Signal - Q1\_Data1')

subplot(2,3,2)
plot(time_Q1_data_2, ERP_Q1_Data2)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q1\_Data2')
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
plot([0 0]+.5,get(gca,'ylim'),'k--')

subplot(2,3,5)
plot(hz_Q1_data_2, abs(fft(ERP_Q1_Data2)),"linew",2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectrum of Average ERP Signal - Q1\_Data2')

subplot(2,3,3)
plot(time_Q1_data_3, ERP_Q1_Data3)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q1\_Data3')
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
plot([0 0]+.5,get(gca,'ylim'),'k--')


subplot(2,3,6)
plot(hz_Q1_data_3, abs(fft(ERP_Q1_Data3)),"linew",2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectrum of Average ERP Signal - Q1\_Data3')


%-------------------------------------------------------Q2.1---------------------------------------

% Loading data file
Q2_Data = load('Q2_Data.mat');
% Extracting first channel
Q2_Data_channel1 = Q2_Data.EEG.data(1, :, :); % first channel for Q2_Data
% Finding average ERP signal using averaging
ERP_Q2_Data = mean(Q2_Data_channel1, 3); % average ERP signal for Q2_Data
% Seting time axis for plot
time_Q2 = 0:1/Q2_Data.EEG.srate:(length(ERP_Q2_Data)-1)/Q2_Data.EEG.srate; % time axis for Q2_Data
% Ploting average ERP signal in time domain
figure
plot(time_Q2, ERP_Q2_Data)
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
plot([0 0]+.5,get(gca,'ylim'),'k--')
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q2\_Data (Channel 1)')

%----------------------------------------------------------Q2.2----------------------------------------


% Loading data file
Q2_Data = load('Q2_Data.mat');

% Extracting first channel
Q2_Data_channel1 = Q2_Data.EEG.data(1, :, :); % first channel for Q2_Data
% Finding average ERP signal in time domain
ERP_Q2_Data = mean(Q2_Data_channel1, 3); % average ERP signal in time domain
% Seting time axis for plot
time_Q2 = 0:1/Q2_Data.EEG.srate:(length(ERP_Q2_Data)-1)/Q2_Data.EEG.srate; % time axis for Q2_Data
% Calculating power spectrum of first trial
trial1_fft_Q2 = fft(Q2_Data_channel1(:,:, 1)); % Fourier Transform of first trial
trial1_power_spectrum_Q2 = abs(trial1_fft_Q2).^2/length(trial1_fft_Q2); % power spectrum of first trial

frequencyRepSeparate_Q2 = fft(squeeze(Q2_Data_channel1(:,:,1)))/length(time_Q2);
powspectSeparate_Q2 = mean((2*abs(frequencyRepSeparate_Q2)).^2,2);
% Calculating power spectrum of averaged trials in time domain
time_domain_average_fft_Q2 = fft(ERP_Q2_Data); % Fourier Transform of average ERP signal in time domain
time_domain_average_power_spectrum_Q2 = abs(time_domain_average_fft_Q2).^2/length(time_domain_average_fft_Q2); % power spectrum of average ERP signal in time domain
% Calculating power spectrum by averaging the Fourier representations of individual trials
trials_fft_Q2 = fft(Q2_Data_channel1(:,:,:)); % Fourier Transform of all trials      , [], 2
trials_power_spectrum_Q2 = mean(abs(trials_fft_Q2).^2, 3)/size(trials_fft_Q2, 2); % power spectrum by averaging the Fourier representations of individual trials
% our hz axis
hz_Q2=(0:1200-1)*400/1200;
% Ploting powear spectrums
figure
subplot(3, 1, 1)
plot(hz_Q2, trial1_power_spectrum_Q2,"linew",2)
xlim([0 100])
title('Power Spectrum of First Trial')

subplot(3, 1, 2)
plot(hz_Q2, time_domain_average_power_spectrum_Q2,"linew",2)
xlim([0 100])
title('Power Spectrum of Averaged Trials in Time Domain')

subplot(3, 1, 3)
plot(hz_Q2, trials_power_spectrum_Q2,"linew",2)
xlim([0 100])
title('Power Spectrum Calculated by Averaging the Fourier Representations of Individual Trials')

% Loading data file
%--------------------------------------------Q3.1----------------------------------------------
Q3_data = load('Q3_Data.mat');
sampling_rate = Q3_data.EEG.srate;
Q3_data_channels = Q3_data.EEG.data(1, :, :);
ERP_Q3_data = mean(Q3_data_channels,3);
tQ3 = 0:1/Q3_data.EEG.srate:(length(ERP_Q3_data)-1)/Q3_data.EEG.srate;

figure
subplot(2,1,1)
plot(tQ3, ERP_Q3_data,"linew",1.25)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q3/data')

[a,b] = xcorr(ERP_Q3_data,ERP_Q3_data); % Autocorrelation for center frequency 
subplot(2,1,2)
plot(b,a,"linew",1.25)
xlabel('Autocorelated f center spect')
ylabel('F value')
title('Centered Freq Spect')

Bandwidth = rms(ERP_Q3_data); % Root mean square for bandwith
Bandwidth

%-----------------------------------------------Q3.2------------------------------------------------

trial1_fft_Q3_data = fft(Q3_data_channels(:,:, 1));
trial1_power_spectrum_Q3_data = abs(trial1_fft_Q3_data).^2/length(trial1_fft_Q3_data);

frequencyRepSeparate_Q3_data = fft(squeeze(Q3_data_channels(:,:,1)))/length(tQ3);
powspectSeparate_Q3_data = mean((2*abs(frequencyRepSeparate_Q3_data)).^2,2);

time_domain_average_fft_Q3_data = fft(ERP_Q3_data);
time_domain_average_power_spectrum_Q3_data = abs(time_domain_average_fft_Q3_data).^2/length(time_domain_average_fft_Q3_data);

trials_fft_Q3_data = fft(Q3_data_channels(:,:,:));
trials_power_spectrum_Q3_data = mean(abs(trials_fft_Q3_data).^2, 3)/size(trials_fft_Q3_data, 2);

hz_Q3_data=(0:1200-1)*400/1200;

figure

subplot(3, 1, 1)
plot(hz_Q3_data, trial1_power_spectrum_Q3_data,"linew",1.25)
xlim([0 100])
title('Power Spectrum of First Trial')

subplot(3, 1, 2)
plot(hz_Q3_data, time_domain_average_power_spectrum_Q3_data,"linew",1.25)
xlim([0 100])
title('Power Spectrum of Averaged Trials in Time Domain',"linew",1.25)

subplot(3, 1, 3)
plot(hz_Q3_data, trials_power_spectrum_Q3_data,"linew",1.25)
xlim([0 100])
title('Power Spectrum Calculated by Averaging the Fourier Representations of First Trials')

%--------------------------------------------Q4.1---------------------------------------------
%loading data
Q4_Data = load('Q4_Data.mat');

% Taking first channel data
Q4_Data_channel1 = Q4_Data.EEG.data(1, :, :);

% ERP for each trials
ERP_Q4_Data = mean(Q4_Data_channel1, 3);

% Making time axis
time_Q4 = 0:1/Q4_Data.EEG.srate:(length(ERP_Q4_Data)-1)/Q4_Data.EEG.srate;


% spectrogram f-graph making
figure
subplot(2,3,1)
spectrogram(ERP_Q4_Data, hamming(256), 240, 1000, 400);
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of Average ERP Signal with Hamming Window')

subplot(2,3,2)
spectrogram(ERP_Q4_Data, hann(256), 240, 1000, 400);
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of Average ERP Signal with Hann Window')


subplot(2,3,3)
spectrogram(ERP_Q4_Data, gausswin(256), 240, 1000, 400);
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of Average ERP Signal with Gausian Window')


subplot(2,3,4)
spectrogram(ERP_Q4_Data, rectwin(256), 250, 1000, 400);
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of Average ERP Signal with Rectangular Window')

subplot(2,3,5)
spectrogram(ERP_Q4_Data, blackman(250), 240, 1000,400);
colorbar;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Spectrogram of Average ERP Signal with Blackman Window')


subplot(2,3,6)
plot(time_Q4, ERP_Q4_Data,"linew",1.25)
xlabel('Time (s)')
ylabel('Amplitude (\muV)')
title('Average ERP Signal in Time Domain - Q4')
hold on 
plot(time_Q4, Q4_Data_channel1(1,:,1))

hz_Q4=(0:1200-1)*400/1200;


[~, time_idx] = max(abs(ERP_Q4_Data)); % Find the index of the time point with the maximum absolute value
time_point_Q4 = time_Q4(time_idx) % Convert the index to the time point in seconds

%finding the peak frequency of the transient activity, you can use the following code:

[~, frequency_idx_Q4] = max(abs(ERP_Q4_Data)); % Find the index of the frequency with the maximum absolute value
peak_frequency_Q4 = hz_Q4(frequency_idx_Q4) % Convert the index to the frequency in Hz

%finding the bandwidth of the transient activity, you can use the following code:

[~, frequency_idx_Q4] = max(abs(ERP_Q4_Data)); % Find the index of the frequency with the maximum absolute value
bandwidth_Q4 = hz_Q4(frequency_idx_Q4) - hz_Q4(1) % Calculate the bandwidth by subtracting the lowest frequency from the peak frequency

