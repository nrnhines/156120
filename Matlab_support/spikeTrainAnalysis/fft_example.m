%file to comput the fourier spectrum of my signal 
%using the pwelch function.
Fs=fs;
win_PSD=5*Fs;
noverlap_PSD=[];
nFFT=2.^10;
[spd,f]=pwelch(Outp,win_PSD,noverlap_PSD,nFFT,Fs); %spd=spectral power
                                                  %arbitrary units/Hz
