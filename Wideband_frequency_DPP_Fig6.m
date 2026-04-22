clear all;
clc;
close all;

%% UCA settings
Nt = 256;
fc = 30e9;
c = 3e8;
lambda = c/fc;
d = lambda/2;

r_radius = Nt*d/2/pi; % half-wavelength spacing
D = 2*(2*r_radius)^2/lambda;
C = linspecer(6);

%% subarray response analysis
% P = 2;
% Nb = Nt/P;
% 
% pp = 1:1:Nt;
% phi_p = (pp-1)*2*pi/Nt;
% theta0 = 0;
% 
% dis_matrix = -r_radius*cos(theta0-phi_p);
% at = exp(-1j*2*pi/lambda*dis_matrix);
% sum_all = sum(at);
% sum_all_est = Nt*abs(besselj(0, 2*pi*r_radius/lambda));
% 
% sum_sub = sum(at(1:Nb));
% M_bessel = 1000;
% M_bessel_list = -M_bessel:1:M_bessel;
% sub_est = zeros(1, length(M_bessel_list));
% Nb_list = 0:Nb-1;
% for i_m = 1:length(M_bessel_list)
%     m = M_bessel_list(i_m);
%     geometric_sum = sum(exp(-1j*Nb_list*m*2*pi/Nt));
%     sub_est(i_m) = ((1j)^m)*besselj(m, 2*pi*r_radius/lambda)*geometric_sum;
% end
% sum_sub_est = sum(sub_est);
% 
% ratio = sum_sub/sum_sub_est;

%% wideband gain estimation
% classical narrow band
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;
phi0 = 0;
dis_matrix = -r_radius*cos(phi0-phi_p);
at = exp(-1j*2*pi/lambda*dis_matrix)/sqrt(Nt);

% bandwidth
BW = 3e9;
num_carrier = 801;
f_bw = linspace(fc-BW/2, fc+BW/2, num_carrier);
at_bw = exp(-1j*2*pi*f_bw.'/c*dis_matrix)/sqrt(Nt);

% No DPP
gain_bw = abs(at_bw*at');
figure;
hold on;
plot(f_bw/1e9, gain_bw,'-','LineWidth',1.5,'color',C(1,:));
xlabel('f');
ylabel('Gain');
grid on;
box on;


% subarray
K = 8;  % number of sub-arrays
P = Nt/K;
A_phase = zeros(Nt, K);
for i_k = 1:K
    A_phase((i_k-1)*P+1:i_k*P, i_k) = at((i_k-1)*P+1:i_k*P);
end
pp_time = 1:1:K;
phi_p_time = (2*(pp_time-1)+1)*pi/K;

dis_matrix_time = -r_radius*cos(phi0-phi_p_time);
at_time = exp(1j*2*pi*(fc-f_bw).'/c*dis_matrix_time);

gain_bw_time = abs(sum(A_phase*at_time.'.*at_bw'));

% %% infinite summation
% Nb_list = 0:Nb-1;
% M_bessel = 8;
% M_bessel_list = -M_bessel:1:M_bessel;
% gain_bw_record = zeros(num_carrier, 2*M_bessel+1);
% for i_bw = 1:num_carrier
%     fk = f_bw(i_bw);
%     for i_m = 1:length(M_bessel_list)
%         m = M_bessel_list(i_m);
%         geometric_sum = sum(exp(-1j*Nb_list*m*2*pi/Nt));
%         gain_bw_record(i_bw, i_m) = 1/Nb*besselj(m, 2*pi*r_radius/c*(fc-fk))^2 ...
%             *exp(1j*m*pi/P)*geometric_sum;
%     end
% end
% gain_bw_est = abs(sum(gain_bw_record,2));
% plot(f_bw, gain_bw_est,'b--','LineWidth',1.5);

%% estimation with one item
% R = sqrt(2)*2*pi*r_radius/c*(fc-f_bw)*sqrt(1-cos(pi/P/2+pi/Nt/2));
% gain_bw_zero = abs(besselj(0, R));
% plot(f_bw, gain_bw_zero,'b--','LineWidth',1.5);

% gain_bw_zero = abs(besselj(0, 2*pi*r_radius/c*(fc-f_bw)).^2);
% plot(f_bw, gain_bw_zero,'k--','LineWidth',1.5);

%% estimation with J0
gain_bw_re = zeros(num_carrier, P);
for i_nb = 1:P
    i = i_nb-1;
    Ri = sqrt(2)*2*pi*r_radius/c*(fc-f_bw)*sqrt(1-cos(pi/K-i*2*pi/Nt));
    gain_bw_re(:, i_nb) = besselj(0, Ri);
end
gain_bw_sum = abs(sum(gain_bw_re, 2))/P;

%% hypergeometric function
a = 1/2;
b = [1,3/2];

gain_bw_sum_hyper = zeros(1, num_carrier);
for i_bw = 1:num_carrier
    fk = f_bw(i_bw);
    zeta = P*pi*r_radius/c*(fc-fk)*2*pi/Nt;
    z = -zeta^2/4;
    gain_bw_sum_hyper(i_bw) = hypergeom(a,b,z);
end

idx = 1:100:num_carrier;
figure;
hold on;
h_1 = plot(f_bw/1e9, gain_bw,'-','LineWidth',1.5,'color',C(1,:));

h2 = plot(f_bw/1e9, gain_bw_sum,'--','LineWidth',1.5,'color',C(3,:));
h2_1 = plot(f_bw(idx)/1e9, gain_bw_sum(idx),'+','LineWidth',1.5,'color',C(3,:));
h2_2 = plot(f_bw(1)/1e9, gain_bw_sum(1),'--+','LineWidth',1.5,'color',C(3,:));

h3 = plot(f_bw/1e9, gain_bw_time,'--','LineWidth',1.5,'color',C(2,:));
h3_1 = plot(f_bw(idx)/1e9, gain_bw_time(idx),'o','LineWidth',1.5,'color',C(2,:));
h3_2 = plot(f_bw(1)/1e9, gain_bw_time(1),'--o','LineWidth',1.5,'color',C(2,:));

h4 = plot(f_bw/1e9, gain_bw_sum_hyper,'--','LineWidth',1.5,'color',C(4,:));
h4_1 = plot(f_bw(idx)/1e9, gain_bw_sum_hyper(idx),'d','LineWidth',1.5,'color',C(4,:));
h4_2 = plot(f_bw(1)/1e9, gain_bw_sum_hyper(1),'--d','LineWidth',1.5,'color',C(4,:));

xlabel('Frequency (GHz)');
ylabel('Gain');
grid on;
box on;
% print_str = ['K=', num2str(K)];
% title(print_str);
h=legend([h_1,h2_2,h3_2,h4_2], 'No DPP','DPP Accurate','Estimation with (14)','Estimation with (15)');
set(h,'Location','northeast');
% legend('No DPP','DPP Accurate','Estimation with (13)','Estimation with (14)');


% idx = 1:10:num_carrier;
% idx = [idx, plot_point1:5:plot_point1+plot_point2];
% h1 = plot(D_list,abs(N_cal),'-','color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
% h1_1 = plot(D_list(idx),abs(N_cal(idx)),'^','color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
% h1_2 = plot(D_list(1),abs(N_cal(1)),'-^','color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
% h2 = plot(D_list,abs(N_est1),'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5);
% h2_1 = plot(D_list(idx),abs(N_est1(idx)),'x','color',[0.4940 0.1840 0.5560],'LineWidth',1.5);
% h2_2 = plot(D_list(1),abs(N_est1(1)),'-x','color',[0.4940 0.1840 0.5560],'LineWidth',1.5);
% h3 = plot(D_list,abs(N_est2),'-','color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
% h3_1 = plot(D_list(idx),abs(N_est2(idx)),'o','color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
% h3_2 = plot(D_list(1),abs(N_est2(1)),'-o','color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
% h=legend([h1_2,h2_2,h3_2], 'Channel DoF', 'Rough Estimation', 'Accurate Estimation', 'FontSize',13);
% set(h,'Interpreter','latex','Location','northeast');




