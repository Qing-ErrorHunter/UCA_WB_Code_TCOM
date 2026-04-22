clear all;
clc;
close all;

% CoreNum=32;
% if isempty(gcp('nocreate'))
%     parpool(CoreNum);
% end

%% UCA settings
N = 256;
fc = 30e9;
c = 3e8;
M=4;
lambda = c/fc;
d = lambda/2;
r_radius = N*d/2/pi;
L = 4;
tmax=20e-9;

BW = 3e9;
num_sub_car = 128;
SNR_dB= 0;
SNR=10.^(SNR_dB/10);

K_list = 2.^(0:1:5);
num_k = length(K_list);

N_iter = 1; % increase the number of iterations to obtain a smooth curve

a = 1/2;
b = [1,3/2];

rate_optimal = zeros(num_k, N_iter);
rate_optimization = zeros(num_k, N_iter);
rate_TTD = zeros(num_k, N_iter);
rate_spatial = zeros(num_k, N_iter);
rate_estimated = zeros(num_k, N_iter);
rate_all = zeros(num_k, N_iter);
rate_multiple  = zeros(num_k, N_iter);
for i_k=1:num_k
    K = K_list(i_k);
    P = N/K;
    for i_iter =1:N_iter
        [H,hc,Theta,Alpha, beta] = Wideband_Channel_UCA_Multi(N,M,L,fc,BW,tmax,num_sub_car);
         
        %% Optimal
        F_max = zeros(num_sub_car,N);
        for j=1:num_sub_car
            [~,~,V]=svd(H(:,:,j));
            F = V(:,1:M);
            F = F/norm(F);
            HH = H(:,:,j)*F;
            rate_optimal(i_k, i_iter) = rate_optimal(i_k, i_iter)...
                +real(log2(det(eye(M)+HH*HH'*SNR)));
        end

        %% Hybrid near optimal
        R=zeros(N,N);
        for j=1:num_sub_car
            R=R+1/num_sub_car*H(:,:,j)'*H(:,:,j);
        end
        [~,~,V]=svd(R);
        F1_temp = V(:,1:M);
        for j=1:num_sub_car
            H1_temp = H(:,:,j)*F1_temp;
            [~,~,Vbb] = svd(H1_temp);
            Fbb = Vbb(:,1:M);
            F1 = F1_temp*Fbb/norm(F1_temp*Fbb);
            HH1 = H(:,:,j)*F1;
            rate_optimization(i_k, i_iter) = rate_optimization(i_k, i_iter)...
                +real(log2(det(eye(M)+HH1*HH1'*SNR)));
        end 

        %% Partial TTD
        [F_PS, t_k_list] = DPP_uca(N,min(L,M),K,r_radius,fc,Theta);
        for j=1:num_sub_car
            f=fc+BW/(num_sub_car)*(j-1-(num_sub_car-1)/2);
            F_TTD = zeros(K*min(L,M),min(L,M));
            for i_m=1:min(L,M)
                F_TTD((i_m-1)*K+1:i_m*K, i_m) = exp(-1j*2*pi*f*t_k_list(:,i_m));
            end
            F3_temp = F_PS*F_TTD;
            H3_temp = H(:,:,j)*F3_temp;
            [~,~,Vbb] = svd(H3_temp);
            Fbb = Vbb(:,1:M);
            F3 = F3_temp*Fbb/norm(F3_temp*Fbb);
            HH3 = H(:,:,j)*F3;
            rate_TTD(i_k, i_iter) = rate_TTD(i_k, i_iter)+...
                real(log2(det(eye(M)+HH3*HH3'*SNR)));
        end
        
        %% estimated
        for i_bw = 1:num_sub_car
            fk=fc+BW/num_sub_car*(i_bw-1-(num_sub_car-1)/2);
            zeta = P*pi*r_radius/c*(fc-fk)*2*pi/N;
            z = -zeta^2/4;
            gain_hyper = hypergeom(a,b,z);
            [~,S,V] = svd(squeeze(H(:,:,i_bw)));
            F = V(:,1:M);
            F = F/norm(F);
            HH = H(:,:,i_bw)*F;
            rate_estimated(i_k, i_iter) = rate_estimated(i_k, i_iter)...
                +real(log2(det(eye(M)+S(1:M,1:M).^2*gain_hyper^2*SNR)));
        end

        %% Spatially sparse precoding
        F=[];
        for i_m=1:M
            at = sqrt(1/N)*exp(1j*2*pi*r_radius/c*fc*cos(Theta(i_m)-([0:N-1]*2*pi/N))).';
            F=[F at];
        end
        for j = 1:num_sub_car
            H4_temp = H(:,:,j)*F;
            [~,~,Vbb] = svd(H4_temp);
            Fbb = Vbb(:,1:M);
            F4 = F*Fbb/norm(F*Fbb);
            HH4 = H(:,:,j)*F4;
            rate_spatial(i_k, i_iter) = rate_spatial(i_k, i_iter)+...
                real(log2(det(eye(M)+HH4*HH4'*SNR)));
        end
        
        %% array vectors of all subcarriers (AV-All)
        F_all=[];
        for i_m=1:M
            at_add = zeros(N, 1);
            for i_sub = 1:num_sub_car
                f=fc+BW/(num_sub_car)*(i_sub-1-(num_sub_car-1)/2);
                at_add = at_add + sqrt(1/N)*exp(1j*2*pi*r_radius/c*f*cos(Theta(i_m)-([0:N-1]*2*pi/N))).';
            end
            F_all = [F_all, at_add];
        end
        for j = 1:num_sub_car
            H4_temp = H(:,:,j)*F_all;
            [~,~,Vbb] = svd(H4_temp);
            Fbb = Vbb(:,1:M);
            F_all_norm = F_all*Fbb/norm(F_all*Fbb);
            HH_all = H(:,:,j)*F_all_norm;
            rate_all(i_k, i_iter) = rate_all(i_k, i_iter)+real(log2(det(eye(M)+HH_all*HH_all'*SNR)));
        end
        
        %% array vectors of multiple subcarriers (AV-Multiple)
        F_all=[];
        sub_car_loc = [num_sub_car/4, num_sub_car/4*3];
        num_loc = length(sub_car_loc);
        for i_m=1:M
            at_add = zeros(N, 1);
            for i_sub = 1:num_loc
                sub_idx = sub_car_loc(i_sub);
                f=fc+BW/(num_sub_car)*(sub_idx-1-(num_sub_car-1)/2);
                at_add = at_add + sqrt(1/N)*exp(1j*2*pi*r_radius/c*f*cos(Theta(i_m)-([0:N-1]*2*pi/N))).';
            end
            F_all = [F_all, at_add];
        end
        for j = 1:num_sub_car
            H4_temp = H(:,:,j)*F_all;
            [~,~,Vbb] = svd(H4_temp);
            Fbb = Vbb(:,1:M);
            F_all_norm = F_all*Fbb/norm(F_all*Fbb);
            HH_all = H(:,:,j)*F_all_norm;
            rate_multiple(i_k, i_iter) = rate_multiple(i_k, i_iter)+real(log2(det(eye(M)+HH_all*HH_all'*SNR)));
        end

        fprintf('i_SNR=%d ,i_loop=%d\n',i_k,i_iter);
    end
end

Rate = sum(rate_optimal,2)/num_sub_car/N_iter;
Rate2 = sum(rate_estimated,2)/num_sub_car/N_iter;
Rate3 = sum(rate_TTD,2)/num_sub_car/N_iter;
Rate1 = sum(rate_optimization,2)/num_sub_car/N_iter;
Rate4 = sum(rate_spatial,2)/num_sub_car/N_iter;
Rate_all = sum(rate_all,2)/num_sub_car/N_iter;
Rate_multiple = sum(rate_multiple,2)/num_sub_car/N_iter;

%% Figure
C = linspecer(9);
figure;
hold on;
plot(K_list,Rate,'--','Linewidth',2.0,'Markersize',5,'color', C(2,:));
plot(K_list,Rate3,'-o','Linewidth',1.5,'Markersize',5,'color', C(1,:));
plot(K_list,Rate1,'->','Linewidth',1.5,'Markersize',5,'color', C(4,:));
plot(K_list,Rate4,'-d','Linewidth',1.5,'Markersize',5,'color', C(9,:));
plot(K_list,Rate_multiple,'-*','Linewidth',1.2,'Markersize',5,'color', C(7,:));

grid on;
box on;
l1=xlabel('Number of TTD units');
l2=ylabel('Spectrum Efficiency (bit/s/Hz)');
l3=legend('Optimal fully-digital precoding', 'Proposed DPP algorithm',...
    'Wideband optimizaiton method [15]','Selecting single steering vector [16]',...
    'Combining multiple steering vectors [16]','Location','Northwest');
set(l1,'FontSize',12);
set(l2,'FontSize',12);
set(l3,'FontSize',10);
