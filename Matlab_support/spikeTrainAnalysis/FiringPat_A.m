clear a b b1

total_time=40000;%(ms)              %=========
time_step=6;%(ms)                   %========= Constants Definition
N_elem=round(total_time/time_step); %=========
t_simul=total_time/1000;            %========= Simulated period in seconds
tstep=t_simul/N_elem;               %========= time step in seconds
time_vect=tstep*(1:N_elem);         %=========Time vector creation
a=zeros(1,N_elem);



%Collapsing the whole cell object of spikes
ii=0;
for r=1:size(Raster,1)
    for count=1:size(Raster{r})ii=ii+1; b(ii)=Raster{r}(count);end
end

%Creating the Histogram
b=sort(b);
for i=time_step:time_step:max(b)
    if(i-time_step<b(1)&&b(1)<=i)i_top=i; break; end
end
ii=1;
while(ii<=size(b,2)&&i_top<max(b))
    elem=i_top/time_step;
    while(b(ii)<i_top)a(elem)=a(elem)+1; ii=ii+1; end
    i_top=i_top+time_step;    
end


%Creating the Histogram2 second time step

for i=time_step:time_step:max(b)
    if(i-time_step<=b(1)&&b(1)<=i)i_val=floor(i); break; end
end
c2=zeros(1,N_elem);
i_val0=i_val/1000;%%first value to make the hist in seconds

elem=round(i_val/time_step); epsil=3; i=1;
while(elem<N_elem&&i_val<max(b))   
   while(ii<size(b,2)&&b(ii)<(i_val-epsil)) ii=ii+1;end 
   while(ii<size(b,2)&&((i_val-epsil)<=b(ii)&&b(ii)<=(i_val+epsil)))c2(elem)=c2(elem)+1; ii=ii+1;end
   i_val=i_val+time_step;ii=1;elem=elem+1;
 end




