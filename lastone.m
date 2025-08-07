cd('/storage2/teaching/practicum_evgeny/data_analysis/users/ephys01/ER02-20110919/')
mfilename = 'Practicum_Theta';
FileBase = 'ER02-20110919';
ThetaChannel = 72;
BrainState = 'RUN';
SignalType = 'lfp';

par = LoadXml('ER02-20110919.xml');
lfpSamplingRate = par.lfpSampleRate;
BhvStatePeriods = dlmread('ER02-20110919.sts.RUN', '\t');
lfp = LoadBinary('ER02-20110919.lfp', ThetaChannel);
lfp = lfp';
spd = load('ER02-20110919.spd');

%plot(spd)

lfpSamplingRate = 1250;
VideoSamplingRate = 50;

BhvStatePeriodsLFP = BhvStatePeriods;
lfp = SelectPeriods(lfp, BhvStatePeriodsLFP, 'c');

BhvStatePeriods_sec = BhvStatePeriodsLFP / lfpSamplingRate;
BhvStatePeriodsSpeed = round(BhvStatePeriods_sec * VideoSamplingRate);
spd = SelectPeriods(spd, BhvStatePeriodsSpeed, 'c');

%DEBUGGING FIGURE
figure;
%Plot LFP for a single channel with time in samples
subplot(3,1,1); 
% plot(lfp, 'Color','k','LineStyle','-', 'LineWidth',0.5 );
t_lfp = [0 : size(lfp,1)-1] / lfpSamplingRate;
plot(t_lfp, lfp, 'Color','k','LineStyle','-', 'LineWidth',0.5 );
axis tight
xlabel('Time, (sec)')
ylabel('Amplitude')
title('LFP signal')
subplot(3,1,2)
t_spd = [0 : length(spd)-1] / VideoSamplingRate;
plot(t_spd, spd, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 0.5);
axis tight
xlabel('Time, (sec)')
ylabel('Speed (cm/s)')
title('Animal speed')

% to this point its revised and fine
% 
% % Duration to zoom in (sec)
% t_zoom = 10;
% 
% % Index range for zoomed LFP and speed
% lfp_range = 1 : round(t_zoom * lfpSamplingRate);
% spd_range = 1 : round(t_zoom * VideoSamplingRate);
% 
% % Time vectors
% t_lfp_zoom = lfp_range / lfpSamplingRate;
% t_spd_zoom = spd_range / VideoSamplingRate;
% 
% % Plot figure
% figure;
% 
% % === LFP ZOOMED ===
% subplot(2,1,1);
% plot(t_lfp_zoom, lfp(lfp_range), 'k-', 'LineWidth', 1);
% xlabel('Time (sec)');
% ylabel('LFP amplitude');
% title('Zoomed LFP (theta visible)');
% axis tight
% ylim([-10000 10000]);  % Adjust based on your signal amplitude to make theta visible
% 
% % === SPEED ZOOMED ===
% subplot(2,1,2);
% plot(t_spd_zoom, spd(spd_range), 'b-', 'LineWidth', 1);
% xlabel('Time (sec)');
% ylabel('Speed (cm/s)');
% title('Corresponding Animal Speed');
% axis tight
% ylim([0 max(spd(spd_range)) + 5]);  % Add padding above peak speed



%Set parameters for multitaper spectral calculations 
%Frequency range the power spectrum must be computed over (Hz)
FreqRange = [1 20];

%Length of moving spectral window (sec)
%It should accomodate at least few cycles of the oscillation of interest)
WinLength_sec = 3;

%Convert the window length from sec to sample indices, 
%because the function which computes PSD requires it.
WinLength = round(WinLength_sec*lfpSamplingRate);


%Overlap between the successive windows, expressed as a fraction of the window length 
%(set to 80-90% for good smoothing over time)
nOverlap = round(0.9*WinLength);

%"Time-bandwidth product" for the taper functions, controling the degree of smoothing of PSD profile.
%(frequency resolution)
%The function mtcsdfast.m uses the first 2*NW-1 tapers for the in the spectral estimate (nTapers = 2*NW-1).
NW = 2;

%Number of points of in the FFT to calculate.
%nFFT controls the number of frequency bins covering the frequency range of interest
%If nFFT is greater than the signal length, the signal is zero-padded to the length nFFT.
%If nFFT is less than the signal length, the signal is cut to the length nFFT.
nFFT = 2000;



%Compute spectrogram (time-frequency decomposition)
[spec.y, spec.f, spec.t] = mtcsglong(lfp, nFFT, lfpSamplingRate, WinLength, nOverlap, NW, [],[], FreqRange);


%OPTIONAL: Swap the dimensions of spec.y now, to avoid doing it while plotting below
spec.y = spec.y';

%OPTIONAL: Convert power into log10 scale (dB/Hz) for better visualization (it squeezes the power values range)
%Replace zero values with the smallest double precision number in order to avoid -Inf vlues after log.
spec.y(spec.y==0) = eps; 
spec.y = 20*log10(spec.y);

%Plot spectrogram (dB/Hz) 
subplot(3,1,3); cla; hold on
colormap(jet) 
% colormap(parula) 
imagesc(spec.t, spec.f, spec.y);
axis tight; axis xy;
set(gca, 'ylim', [1 20], 'ytick', 0:5:20)
xlabel('Time, (s)')
ylabel('Frequency, (Hz)')
colorbar
title('Spectrogram, (dB/Hz)')



%--------------------------------------------------------------------------------------------%
%              Quantification of theta frequency and power
%                      based on the spectrogram
%--------------------------------------------------------------------------------------------%
%Set necessary parameters
%Frequency range of theta oscillation (Hz)
FreqRangeTheta = [4 10];
%Frequency range for power normalization (Hz)
FreqRangeNorm = [1 20];

       
%Indices of spectrogram frequency bins from the theta range 
SpecFreqBinsTheta = find( spec.f >= FreqRangeTheta(1) & spec.f <= FreqRangeTheta(2) );

%Indices of spectrogram frequency bins from the normalization range
SpecFreqBinsNorm = find( spec.f >= FreqRangeNorm(1) & spec.f <= FreqRangeNorm(2) );

%Extract and sum up the spectrogram power values from the theta frequency range 
ThetaPowerRaw = mean( spec.y(SpecFreqBinsTheta,:), 1)';

%Extract and sum up the spectrogram power values from the normalization frequency range 
NormPowerRaw  = mean( spec.y(SpecFreqBinsNorm,:), 1)';

%Normalize power from the theta range by the power from the normalization range
raw.ThetaPower = ThetaPowerRaw ./ NormPowerRaw;


%Compute theta frequency within each time bin, defined as a frequency of the peak power in each time bin
[~, ind] = max( spec.y(SpecFreqBinsTheta,:), [], 1 );
raw.ThetaFreq = spec.f(SpecFreqBinsTheta(ind));
clear ind



%Plot spectrogram + the labeled peak power bins
subplot(3,1,3); cla; hold on
colormap(jet) 
% colormap(parula) 
imagesc(spec.t, spec.f, spec.y);
axis tight; axis xy;
set(gca, 'ylim', [1 20], 'ytick', 0:5:20)
xlabel('Time, (s)')
ylabel('Frequency, (Hz)')
colorbar
title('Spectrogram, (dB/Hz)')
%label the theta freq range
plot(xlim, FreqRangeTheta(1)*[1 1], 'Color','k','LineStyle','--','LineWidth', 1);
plot(xlim, FreqRangeTheta(2)*[1 1], 'Color','k','LineStyle','--','LineWidth', 1);
%label detected peaks in power within the theta range
plot(spec.t, raw.ThetaFreq, 'Color','r','LineStyle','none','Marker', '.');

%------------------next fig----------------------

%Extract speed values for the timestamps spec.t
spec.t_spd = round( spec.t * VideoSamplingRate );
% raw.Speed = spd(spec.t_spd);

%NEW
raw.Speed = spd(spec.t_spd(2:end));
raw.ThetaPower(1) = [];
raw.ThetaFreq(1)  = [];

%Add random jitter to the ThetaFrequencz values for better appearance
raw.ThetaFreq = raw.ThetaFreq + (rand(size(raw.ThetaFreq))-0.5) *4;


%FIGURE
figure
subplot(3,4,1)
cla; hold on
plot(raw.Speed, raw.ThetaPower, 'Color', 0.7*[1 1 1],'LineStyle','none','Marker', '.' )
xlabel('Animal speed, (cm/sec)')
ylabel('Theta power, (norm.)')


subplot(3,4,2)
cla; hold on
plot(raw.Speed, raw.ThetaFreq, 'Color',0.7*[1 1 1],'LineStyle','none','Marker', '.' )
xlabel('Animal speed, (cm/sec)')
ylabel('Theta frequency, (Hz)')
set(gca, 'ylim', [1 15], 'ytick', 0:5:15)


%Correlation (if pval is less than 0.05, then stat. significant)
[rho_pow, pval_pow] = corr(raw.Speed(:), raw.ThetaPower(:));
[rho_freq, pval_freq] = corr(raw.Speed(:), raw.ThetaFreq(:));


%Linear fit
%Goal to find A and B for ThetaPower = A*Speed + B
p_pow = polyfit(raw.Speed, raw.ThetaPower, 1);

%Plot the fitted lines
x_fit = [0 100];
y_fit = polyval(p_pow, x_fit);


subplot(3,4,1)
cla; hold on
plot(raw.Speed, raw.ThetaPower, 'Color', 0.7*[1 1 1],'LineStyle','none','Marker', '.' )
xlabel('Animal speed, (cm/sec)')
ylabel('Theta power, (norm.)')
plot(x_fit, y_fit, 'Color', 'r','LineStyle','-','Marker', 'none', 'LineWidth', 1.5)
set(gca, 'ylim', [0.99 1.1], 'ytick', [0.99  1.1])
title(sprintf('LInear fit: A=%f1.2, B=%f1.2, CorrCoeff = %f1.2', p_pow, rho_pow ))

%It seems that the strength of correlation between speed and theta frequency 
%(not theta power) for the high speed range (50-100 cm/sec) is weaker than  that for the lower speed range 
 %   (0-50 cm/sec). One could compute CorrCoefficient and the linear fit independently for these two speed ranges (clouds).
  %  In this case, you would need to select first the data points from the lower speed range ind = raw.Speed<50;  
   % corr(raw.Speed(ind), raw.ThetaFreq(ind)); and so on.   Then, same for the higher speed range:  ind = raw.Speed>50;  
    %corr(raw.Speed(ind), raw.ThetaFreq(ind)); ....
ind_low = raw.Speed<50;
[rho_freq_low, pval_freq_low] = corr(raw.Speed(ind_low), raw.ThetaFreq(ind_low));
p_low = polyfit(raw.Speed(ind_low), raw.ThetaFreq(ind_low), 1);

ind_high = raw.Speed > 50;
[rho_freq_high, pval_freq_high] = corr(raw.Speed(ind_high), raw.ThetaFreq(ind_high));
p_high = polyfit(raw.Speed(ind_high), raw.ThetaFreq(ind_high), 1);

figure
plot(raw.Speed, raw.ThetaFreq, '.', 'Color', [0.5 0.5 0.5])
xlabel('Speed (cm/s)')
ylabel('Theta frequency (Hz)')
title('Theta frequency vs Speed')


% Plot for low speed (050 cm/s)
subplot(1,2,1); hold on
ind_low = raw.Speed < 50;
plot(raw.Speed(ind_low), raw.ThetaFreq(ind_low), '.', 'Color', [0.5 0.5 1])
p_low = polyfit(raw.Speed(ind_low), raw.ThetaFreq(ind_low), 1);
y_low = polyval(p_low, x_fit);
plot(x_fit, y_low, 'c', 'LineWidth', 2)
xlabel('Speed (cm/s)')
ylabel('Theta frequency (Hz)')
title(sprintf('LOW speed (<50 (cm/s))\nr=%.2f, p=%.3f', rho_freq_low, pval_freq_low))

% Plot for high speed (>50 cm/s)
subplot(1,2,2); hold on
ind_high = raw.Speed > 50;
plot(raw.Speed(ind_high), raw.ThetaFreq(ind_high), '.', 'Color', [1 0.5 0.5])
p_high = polyfit(raw.Speed(ind_high), raw.ThetaFreq(ind_high), 1);
y_high = polyval(p_high, x_fit);
plot(x_fit, y_high, 'm', 'LineWidth', 2)
xlabel('Speed (cm/s)')
ylabel('Theta frequency (Hz)')
title(sprintf('HIGH speed (>50 (cm/s))\nr=%.2f, p=%.3f', rho_freq_high, pval_freq_high))

%Name of the output figure file (without an extention yet)
FileOut = 'ER02-20110919.ewa33';


%Save figure into a JPG file (function from SirotaLab toolbox)
export_fig(FileOut, '-jpeg', '-nocrop', '-transparent', '-m2')

%Save figure into a PDF file (function from SirotaLab toolbox)
export_fig(FileOut, '-pdf', '-nocrop', '-m2')