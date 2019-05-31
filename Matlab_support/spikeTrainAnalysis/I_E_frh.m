%this file only computes the firing rate histogram for each population

for raster=2:-1:1
    if(raster==1)
        Raster=RasterFS0;       
    elseif(raster==2)
        Raster=Raster_P0;              
    end
    
    FiringPat_A    %This line calls FiringPat

   if(raster==1)        a1=a; b1=c2; b0=b/1000; %%i cells
   elseif(raster==2)    a2=a; b2=c2; b00=b/1000; %%e cells
   end   
   clear a c2
 end 
ta2=time_vect; ta1=ta2;
clear time_vect