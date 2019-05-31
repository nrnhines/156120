function [W,t,frq] = Wavelet_1ch(Data,Fs,Frq_low,Frq_high,Frq_step);
%
% Modified klaus.linkenkaer@cncr.vu.nl, 060328.
%
%
%**************************************************************************
% Purpose:
%
% Compute the wavelet transform in the frequency range from 'Frq_low' to
% 'Frq_high' in steps of 'Frq_step'.
%
%**************************************************************************
% Example...
%
%
%**************************************************************************
% Input...
%
% Data     	: data vector.
% Fs	  	: sampling frequency of the INPUT data matrix M.
% Frq_low 	: lowest frequency in the time-frequency (TF) matrix.
% Frq_high 	: highest frequency in the time-frequency matrix.
% Frq_step 	: spacing of frequency lines in TF matrix or TF plot.
%
%
%**************************************************************************
% Output...
%
% W         : the matrix containing the wavelet transform (complex!).
% t         : vector with time points in units of seconds (for plotting).
% frq       : vector with frequencies in units of Hz (for plotting).
%
%**************************************************************************
% Wavelet analysis...

[W,p,s] = wavelet33(Data,1/Fs,1,Frq_step,Frq_low*1.033,Frq_high*1.033);

t = 1/Fs:1/Fs:length(Data)/Fs;
frq = Frq_low:Frq_step:Frq_high;

