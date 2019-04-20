fs = 44100;

%% generate impulse response coefficients from csv

T = 4; % length of impulse response in seconds
decayRate = 1/T; % ranges [0, 1]

b = zeros(1, fs*T);
test1 = load('test1.csv');
times = test1(find(test1 < T));
periodOfFilter = max(times);
b(floor(times*44100)) = 1-periodOfFilter*decayRate*times;

%% generate logarithmic sinesweep and filter
f_0 = 1; % [Hz]
f_1 = 1e4; % [Hz]
T_p = 2; % [s] Duration of sweep
k = (f_1/f_0).^(1/T_p);
t = (0:(fs*T_p-1))/fs;
phase_init = 0;
sinesweep = [sin(phase_init + 2*pi*f_0*(((k.^t)-1)/log(k))) zeros(1, fs*T_p)];

% impulse = [1, zeros(1, fs-1)];
y = filter(b, 1, sinesweep)/nnz(b);

%% load real audio and filter

filename = 'recording.wav';
recording = audioread(filename);

input = [recording; zeros(length(recording), 1)]';
y = filter(b, 1, input)/max(abs(y));%/nnz(b);




%% generate stft frames from impulse responses
% take 2048 point fft
   % use 1025 points: DC and all positive freq
   % centered hanning window
   % hop length of 0
   % build stft from .wav sampled at 22.05kHz
   % need total of 90,000 frames
   % --> determine how many frames are needed for 
   %        one impulse response
   
   % >>> frames.shape
   % (90663, 2049)
   %    each spectrogram is either length 47, 141, 643, 1929
  
% load 'all_frames.npy' from 'matlab_frames 4.mat'

% adds 'frames' variable
%load('matlab_frames 4.mat');


%% plot spectorgram of canne data

% load npy array
load('matlab_frames 4.mat');
spectrogram = frames';


fft_length = 2*1024;
window_length = fft_length;
overlap = 0;

%[~, freq_vec] = fft_plus(spect_in_data_avg(1:window_length)', fs, fft_length);
freq_vec = (-fft_length/2+1:fft_length/2)/fft_length*fs;

time_vector = linspace(1, length(spectrogram'),  floor(length(spectrogram')/(window_length-overlap)));
%freq_vec = linspace(1, fs/2, fft_length/2);
%freq_vec(fft_length/2+1:fft_length*3/4);
figure(103)
image(time_vector, freq_vec(fft_length/2:fft_length*3/4), 20*log10(abs(spectrogram(1:fft_length/2,:))), 'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('time [s]')
ylabel('frequency [Hz]')
title('Spectrogram')
colorbar


%% plot spectrogram of filtered data
%spect_in_data = frames;
%spect_in_data_avg = sum([spect_in_data(1:fs*4) spect_in_data(fs*4+1:fs*8) spect_in_data(fs*8+1:fs*12) spect_in_data(fs*12+1:fs*16) spect_in_data(fs*16+1:fs*20)], 2)/5;
spect_in_data_avg = y';
%spect_in_data_avg = frames';
%spectrogram = frames';

fft_length = 2*1024;
window_length =  fft_length;%fft_length/2;
overlap = 0;%128;
spectrogram = spectrogram_plus(spect_in_data_avg, fs, fft_length, window_length, overlap);

[~, freq_vec] = fft_plus(spect_in_data_avg(1:window_length)', fs, fft_length);

time_vector = linspace(1, length(spect_in_data_avg'),  floor(length(spect_in_data_avg')/(window_length-overlap)))/fs;
%freq_vec = linspace(1, fs/2, fft_length/2);
%freq_vec(fft_length/2+1:fft_length*3/4);
figure(103)
image(time_vector, freq_vec(fft_length/2:fft_length*3/4), 20*log10(abs(spectrogram(1:fft_length/2,:))), 'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('time [s]')
ylabel('frequency [Hz]')
title('Spectrogram')
colorbar
