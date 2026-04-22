clear all;
clc;
close all;

%% UCA settings
Nt = 256;
fc = 30e9;
c = 3e8;
lambda = c/fc;
d = lambda/2;

BW = 3e9;
num_carrier = 2001;
phi0 = pi/6;
Nt_ula = 256;

r_radius = Nt*d/2/pi; % half-wavelength spacing
D = 2*(2*r_radius)^2/lambda;

C = linspecer(6);

%% classical beam split for ULA
% in angular domain
at_ula = array_respones(phi0,Nt_ula,d,lambda);
f_bw = linspace(fc-BW/2, fc+BW/2, num_carrier);
lambda_bw = c./f_bw;
at_ula_bw = array_respones(phi0,Nt_ula,d,lambda_bw);
gain_ula = abs(at_ula_bw'*at_ula);

num_carrier_plot = 5;
f_bw_plot = linspace(fc-BW/2, fc+BW/2, num_carrier_plot);
num_phi = 1001;
phi_list = linspace(asin(0.45), asin(0.55), num_phi);

gain_f_ula = zeros(num_carrier_plot, num_phi);
for i_f = 1:num_carrier_plot
    f_iter = f_bw_plot(i_f);
    lambda_f = c/f_iter;
    for i_phi = 1:num_phi
        phi_iter = phi_list(i_phi);
        at_f_ula = array_respones(phi_iter,Nt_ula,d,lambda_f);
        gain_f_ula(i_f, i_phi) = abs(at_ula'*at_f_ula);
    end
end

figure;
hold on;
for i_plot = 1:(num_carrier_plot+1)/2
    plot(sin(phi_list), gain_f_ula(i_plot,:),'LineWidth',1.5,'color', C(i_plot,:));
end
xlabel('sin(\phi)');
ylabel('Gain');
legend('f_L','f_M','f_c');
grid on;
box on;
axis([0.45 0.55 0 1]);

%% beam defocus for UCA
% in angular domain
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;

dis_matrix = -r_radius*cos(phi0-phi_p);
at = exp(-1j*2*pi/lambda*dis_matrix)/sqrt(Nt);

num_carrier_plot = 5;
f_bw_plot = linspace(fc-BW/2, fc+BW/2, num_carrier_plot);
num_phi = 1001;
phi_list = linspace(asin(0.45), asin(0.55),num_phi);

gain_f = zeros(num_carrier_plot, num_phi);
gain_est = zeros(num_carrier_plot, num_phi);
for i_f = 1:num_carrier_plot
    f_iter = f_bw_plot(i_f);
    for i_phi = 1:num_phi
        phi_iter = phi_list(i_phi);
        dis_matrix_iter = -r_radius*cos(phi_iter-phi_p);
        at_f = exp(-1j*2*pi*f_iter/c*dis_matrix_iter)/sqrt(Nt);
        gain_f(i_f, i_phi) = abs(at*at_f');
    end
    % estimation
    R = 2*pi*r_radius/c*sqrt(f_iter^2+fc^2-2*f_iter*fc*cos(phi_list-phi0));
    gain_est(i_f,:) = abs(besselj(0, R));
end

figure;
hold on;
for i_plot = 1:(num_carrier_plot+1)/2
    plot(sin(phi_list), gain_f(i_plot,:),'LineWidth',1.5,'color', C(i_plot,:));
end
plot(sin(phi_list), gain_est(1:(num_carrier_plot+1)/2,:),'--k','LineWidth',1.5);

xlabel('sin(\phi)');
ylabel('Gain');
legend('f_L','f_M','f_c','Est');
grid on;
box on;
axis([0.45 0.55 0 1]);