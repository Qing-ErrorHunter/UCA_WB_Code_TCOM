clear all;
clc;
close all;

%% parameters
Nt = 256;
fc=30e9; % Frequencey 
c = 3e8;
lambda=c/fc; % wavelegenth;
d=lambda/2;
r_radius = Nt*d/2/pi; % half-wavelength spacing

%% classical wideband beamforming gain -- beam split
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;
[x_p, y_p] = pol2cart(phi_p, r_radius);

phi0 = 0;
dis_matrix = -r_radius*cos(phi0-phi_p);
at = exp(-1j*2*pi/lambda*dis_matrix)/sqrt(Nt);

% beamforming gain for different frequencies at fixed angle
BW = 3e9;
num_carrier = 8001;
f_bw = linspace(fc-BW/2, fc+BW/2, num_carrier);
at_bw = exp(-1j*2*pi*f_bw.'/c*dis_matrix)/sqrt(Nt);
gain_bw = abs(at_bw*at');

figure;
plot(f_bw/1e9, gain_bw,'b-','LineWidth',1.5);
xlabel('Frequency (GHz)');
ylabel('Beamforming Gain');
grid on;
box on;