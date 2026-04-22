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

num_bw = 20;
BW_list = linspace(1e8,5e9,num_bw);

ave_bw = zeros(1, num_bw);
ave_bw_low = zeros(1, num_bw);
ave_bw_up = zeros(1, num_bw);
ave_bw_up_2 = zeros(1, num_bw);

K_list = [4,8,16];
num_k = length(K_list);
ave_dpp = zeros(num_k, num_bw);
ave_dpp_est = zeros(num_k, num_bw);

a = 1/2;
b = [1,3/2];

a2 = [1/2,1/2];
b2 = [1,3/2,3/2];

for i_bw = 1:num_bw

    %% wideband gain estimation
    % classical narrow band
    pp = 1:1:Nt;
    phi_p = (pp-1)*2*pi/Nt;
    phi0 = 0;
    dis_matrix = -r_radius*cos(phi0-phi_p);
    at = exp(-1j*2*pi/lambda*dis_matrix)/sqrt(Nt);

    % bandwidth
    BW = BW_list(i_bw);
    num_carrier = 800;
    f_bw = linspace(fc-BW/2, fc+BW/2, num_carrier);
    at_bw = exp(-1j*2*pi*f_bw.'/c*dis_matrix)/sqrt(Nt);

    % No DPP
    gain_bw = abs(at_bw*at');
    at_bw_est = abs(besselj(0, (2*pi*r_radius)/c*(f_bw-fc)));
    
    % calculate averaged spectral efficiency
    ave_bw(i_bw) = trapz(f_bw, gain_bw)/BW;
    
    z2 = -(BW*pi*r_radius/c)^2;
    ave_bw_up_2(i_bw)  = sqrt(hypergeom(a2,b2,z2));
    ave_bw_up(i_bw) = sqrt(8*c/pi/pi/BW/r_radius);
%     ave_bw_up_cal = trapz(f_bw, sqrt(c/pi/pi/r_radius./abs(f_bw-fc)))/BW;

    z_low = -(pi*r_radius*BW)^2/4/c/c;
    ave_bw_low(i_bw) = hypergeom(a,b,z_low);
%     ave_bw_low_cal = trapz(f_bw, besselj(0, (2*pi*r_radius)/c*(f_bw-fc)))/BW;
    
    %% With DPP 
    for i_k = 1:num_k
        K = K_list(i_k);  % number of sub-arrays
        P = Nt/K;

        %% hypergeometric function

        zeta = P*pi*r_radius/c*(f_bw-fc)*2*pi/Nt;
        z = -zeta.^2/4;
        gain_bw_sum_hyper = hypergeom(a,b,z);

        a2 = [1/2,1/2];
        b2 = [1,3/2,3/2];
        z2 = -(BW*pi^2*r_radius/2/c/K)^2;

        ave_dpp(i_k, i_bw) = trapz(f_bw, gain_bw_sum_hyper)/BW;
        ave_dpp_est(i_k, i_bw) = hypergeom(a2,b2,z2);
    end
end
C = linspecer(6);
BW_list_plot = BW_list/1e9;
figure;
hold on;
l1 = plot(BW_list_plot, ave_bw,'o-','LineWidth',1.5,'color', C(1,:));
% plot(BW_list, ave_bw_up,'--','LineWidth',1.5,'color', C(2,:));
l2 =plot(BW_list_plot, ave_bw_up_2,'--','LineWidth',1.5,'color', C(2,:));
l3 =plot(BW_list_plot, ave_bw_low,'--','LineWidth',1.5,'color', C(3,:));
l4_1 =plot(BW_list_plot, ave_dpp(1,:),'^-','LineWidth',1.5,'color', C(4,:));
l4_2 =plot(BW_list_plot, ave_dpp(2,:),'^-','LineWidth',1.5,'color', C(4,:));
l4_3 =plot(BW_list_plot, ave_dpp(3,:),'^-','LineWidth',1.5,'color', C(4,:));
l5_1 =plot(BW_list_plot, ave_dpp_est(1,:),'k--','LineWidth',1.5);
l5_2 =plot(BW_list_plot, ave_dpp_est(2,:),'k--','LineWidth',1.5);
l5_3 =plot(BW_list_plot, ave_dpp_est(3,:),'k--','LineWidth',1.5);

xlabel('Bandwidth (GHz)');
ylabel('Averaged Beamforming Gain');
grid on;
box on;
print_str = ['K=', num2str(K)];
legend([l1,l2,l3,l4_1,l5_1], 'Averaged gain w/o DPP','Upper bound w/o DPP in (23)',...
    'Lower bound w/o DPP in (24)',...
    'Averaged gain w/ DPP','Estimated averaged gain w/ DPP in (25)');

axis([0 5 0 1.4]);
