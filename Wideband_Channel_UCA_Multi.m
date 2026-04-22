function [H,hc,Theta,Alpha,Beta] = Wideband_Channel_UCA_Multi(N,M,L,fc,fs,tmax,num_sub_car)
%generate the wideband channel with ray-based channel with UCA, M received antenna
%model, coded by Zidong Wu 2023/4/21
%Transmit antenna number N
%Received antenna number M
%path number L
%central frequency fc
%band frequency fs
%subcarrier number num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8;
lambda = c/fc;
d = lambda/2;
r_radius = N*d/2/pi;
H=zeros(M,N,num_sub_car+1);
Theta = rand(L,1)*pi-pi/2;
Alpha = rand(L,1)*pi-pi/2;
Beta = randn(L,1)+randn(L,1)*1j;
% Beta = ones(L,1);

Delay = rand(L,1)*tmax;
 for k=1:num_sub_car+1
     if k==num_sub_car+1
        f=fc;
     else
        f=fc+fs/(num_sub_car)*(k-1-(num_sub_car-1)/2);
     end
     for j=1:L        
       H(:,:,k)=H(:,:,k)+Beta(j)*exp(-1j*2*pi*Delay(j)*f)*array_respones(Alpha(j),M,d,c/f)...
           *array_respones_uca(Theta(j),N,r_radius,c/f)';
     end        
 end
H = H.*sqrt(N*M/L);
hc = H(:,:,num_sub_car+1);
H = H(:,:,1:num_sub_car);
[~,I] = sort(abs(Beta),'descend');
Theta = Theta(I);
Beta = Beta(I);
Alpha = Alpha(I);

end

