function a=array_respones_uca(azimuth,N,r_radius,lambda)
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1j*2*pi*r_radius*cos(azimuth(i)-[0:N-1]*2*pi/N)./lambda.')).'];
end