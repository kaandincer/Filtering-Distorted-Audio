%% function plotPowerSpectrum(y, Fs)
% Scale signal y to have values between -1.0 and 1.0, then plot its power
% spectrum in a new figure window.
%
% y - time sequence
% Fs - sampling rate for y (Hz)
%
% Matthew Lew 10/25/2018

function plotPowerSpectrum(y,Fs)
    [p,f] = periodogram(y/max(abs(y)),hamming(length(y)),[],Fs);
%     [p,f] = pmtm(y/max(abs(y)),[],[],Fs);
    
    figure;
    semilogx(f,20*log10(abs(p)));
    xlim([20 20E3]);        % limit graph to normal human hearing range
%     axis tight
    ylabel('Power/frequency (dB/Hz)')
    xlabel('Frequency (Hz)')
    title('Power Spectrum')
end