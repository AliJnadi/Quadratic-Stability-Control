function F=input_disturbance(t,p,enabled)
F=zeros(p,1);
if enabled&&t>=2&&t<=12
 env=sin(pi*(t-2)/10)^2;
 F(1)=0.1*env*sin(2*pi*0.5*(t-2));F(p/2+1)=0.1*env*cos(2*pi*0.5*(t-2));
end
end
