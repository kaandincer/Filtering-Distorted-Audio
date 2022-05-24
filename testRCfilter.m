%% Case study 3: Circuits as filters
% *ESE 105* 
%
% *Name: FILL IN HERE*

clear;
close all;

%% Part 4: Transfer functions
h = 1/44100;
f = logspace(log10(20),log10(20E3),101);
H_bandpass = zeros(size(f));

for a = 1:length(f)
    % generate inputs
    t = 0:h:2/f(a);
    V_in = 5 * sin(2*pi*f(a)*t);    % Volts
    
    % compute response
    V_bandpass = RCfilter(V_in,h);
    
    % compute transfer function at frequency f
	H_bandpass(a) = max(abs(V_bandpass))/max(abs(V_in));
end

figure;
loglog(f,H_bandpass); box off;
xlim([20 20E3]);    % limit plot to normal human hearing range
ylabel('|H(f)|');
xlabel('f (Hz)');
snapnow

%% Filter a noisy signal

load('handel.mat');
% load('noisyhandel.mat');
% load('apollo11-main-landing.mat');
% load('noisy-apollo11-main-landing.mat');

% set sampling interval to match sampling rate of the audio signal
h = 1/Fs;

% compute signal output from circuit
VsoundFiltered = RCfilter(Vsound,h);

% compare power spectra
plotPowerSpectrum(Vsound,Fs);
plotPowerSpectrum(VsoundFiltered,Fs);

% play original sound
playSound(Vsound,Fs);

% play sound after circuit filter
playSound(VsoundFiltered,Fs);
