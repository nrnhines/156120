% script to write a wavelet from my 
% convoluted alpha function network activity 
if(idx_1==1)Outp=b1;end %%Icells
if(idx_1==0)Outp=b2;end %%Ecells
f_low=0.01;  % lowest frequency
f_high=70;  % highest frequency
f_step=0.1; % frequency step
[W,t,fq]=Wavelet_1ch(Outp,fs,f_low,f_high,f_step);
in_v=Outp;