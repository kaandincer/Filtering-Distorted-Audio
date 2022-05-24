%% function playSound(y, Fs)
% scales the audio signal in y to range from -1.0 to 1.0, then plays the
% entire sequence. Control is returned once the sequence is finished
% playing. 
%
% y - audio time sequence 
% Fs - sampling rate for y (Hz)
%
% Matthew Lew 10/25/2018

function playSound(y, Fs)
obj = audioplayer(y/max(abs(y)), Fs);
playblocking(obj);
end