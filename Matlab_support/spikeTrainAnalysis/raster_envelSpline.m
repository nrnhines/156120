f1=figure()
if(idx_1==1)Raster=RasterFS0; s='I cell index';end
if(idx_1==0)Raster=Raster_P0; s='E cell index';end
subplot(3,1,1);
for i=1:size(Raster,1)
for spike=1:size(Raster{i},1)
x=Raster{i}(spike);
%%plot([x,x+0.3]/1000,[i,i+1.1],'Color','k','LineWidth',0.2); %%here we plot the result against
plot([x,x]/1000,[i,i],'.','MarkerSize',4,'Color','k'); %%here we plot the result against
%%cell member by using little lines; time (sec)
end
set(gca,'FontSize',13);
ylabel(s,'FontSize',14);
xlim([5 15]); %horizontal range
ylim([-2 size(Raster,1)+5])% here we set the highest limit for the rasterplot
hold on;
end

%%Removing exploting splines components due to border conditions
%%=========================================================================
%Splines analyisis
%%=========================================================================


w_comp=200; %wavelet component
%Pretreatment of a data set to bild cubic splines input vector
clear in_v1 thrs0 t thrs1 n_l
t_simul=40; %%total simulated period in seconds
tstep=t_simul/(size(in_v,2));%considering a maximun duration of 40sec

%=======================================================================
i=0;c1=0;
n_l=mean(in_v);%standard noise level in the original signal
clear spk_v0 spk_t0 r
while(i<size(in_v,2)-1)
    i=i+1;        
    if(in_v(i)>n_l)        
        r=1;
        clear a_v
        a_v(r)=in_v(i);        
        while(in_v(i+r)>n_l&&i+r<size(in_v,2))
            r=r+1;
            a_v(r)=in_v(i+r-1);            
        end
        [mx,idx]=max(a_v);
        c1=c1+1;
        spk_v0(c1)=mx;
        spk_t0(c1)=i+idx-1;
        if(r>1)i=i+r;end        
    end
end



%=======================================================================
%this part determines the threshold time_lag assumed periodicity of the
%signal
c1=0;
tOld=0;
for(i=1:size(spk_t0,2))    
        c1=c1+1;
        deltat_V(c1)=spk_t0(i)-tOld;
        tOld=spk_t0(i);
    
end

clear bin_n bin_v h_v
bin_n=200;                          %bins number
bstep=max(deltat_V)/bin_n;          %step between bins
bin_v=bstep/2:bstep:max(deltat_V);  %binned vector
h_v=hist(deltat_V,bin_n)/size(deltat_V,2);         
%thrs1=8
thrs1=6%floor(mean(deltat_V))%%period of the original signal if there is any 
;
%=======================================================================
%This part select the participating elements for the construction of the
%splines coefficients
clear deltat spk_v spk_t 
jitter=0.5;
c1=0;
t=1;
%while(in_v(t)==0)t=t+1;end
while(t<=size(in_v,2))    
    if(in_v(t)>0)        
       clear a_v;
       r1=0; r0=t-floor(jitter*thrs1);
       if(r0<=0)r0=t;end
       for r=r0:t+ceil(jitter*thrs1);
            if(r<size(in_v,2))
                r1=r1+1;                              
                a_v(r1)=in_v(r); 
            end
        end
        [mx,idx]=max(a_v);
        c1=c1+1;
        spk_v(c1)=mx;
        spk_t(c1)=tstep*(t-floor(jitter*thrs1)+idx-1);
        t=t-floor(jitter*thrs1)+idx+thrs1-1;
    else
        t=t+1;
    end
        
end

t_v1=tstep:tstep:t_simul;
splines_V=spline(spk_t,spk_v);
y_val=ppval(splines_V,t_v1);

%%=========================================================================
%Splines analyisis
%%=========================================================================

t=21;
clear tpl1 ttpl1
while(t<size(y_val,2))
tpl1(t-20)=y_val(t);
t=t+1;
end
tpl1=tpl1/1;%%max(tpl1);
ttpl1=21:size(y_val,2)-1;
ttpl1=tstep*ttpl1;


%%=========================================================================


spl_vect=tpl1/1;%%max(tpl1);%%max(in_v);           %%the spline interpolated vector
spl_vect1=abs(W(w_comp,:).^2/max(W(w_comp,:).^2));
spl_vect2=in_v/1;%%max(in_v);
%Y_lim=1.2*max(spl_vect2);
Y_lim=1.2*size(Raster,1);


subplot(3,1,2);
plot(ttpl1,spl_vect,'r','LineWidth',2)%interpolated spline
%hold on
%plot(t_v1,spl_vect1,'LineWidth',2)%200th wavelet component (20Hz)
hold on
bar(t_v1,spl_vect2,'c')
set(gca,'FontSize',13);
%xlim([0 t_simul+0.1]); ylim([0 Y_lim]);
xlim([5 15]); ylim([0 Y_lim]);
ylabel('Firing Cells / 6ms','FontSize',14); 
%%xlabel('time (sec)','FontSize',14);
hold on 
subplot(3,1,3);
[W,t,fq]=Wavelet_1ch(Outp,fs,f_low,f_high,f_step);
imagesc(t,fq,flipud(abs(W)));
colormap(jet(64))
xlim([5 15]); ylim([10 30]);
colorbar('delete')

axis xy;
%%set(gca,'YDir','reverse','FontSize',13)
set(gca,'FontSize',13);
colorbar;
%xlim([5 15]); ylim([0 30]);
xlabel('time (s)')
ylabel('Frequency (Hz)')
set(f1,'Position',[347   376   929   538]);