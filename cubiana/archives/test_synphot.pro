w=ramp(100,5.,12.,/pow)
Fnu=w^2
Fnu+=exp(-0.5*(w-8.)^2/0.1^2)*30          
plotaxis,w,Fnu,psym=-4
Fnu0=synthetic_photometry(w,Fnu,filters=["IRAC3","IRAC4"],wcen=w0)    
opl,w0,fnu0,line=-1,psym="circ",color=!red

end
