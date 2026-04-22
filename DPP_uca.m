function [F_PS, t_k_list] = DPP_uca(N,M,K,r_radius,fc,azimuth)

P = N/K;
c = 3e8;
F_PS = [];
t_k_list = zeros(K,M);
for i_m=1:M
    ac= sqrt(1/N)*exp(1j*2*pi*r_radius/c*fc*cos(azimuth(i_m)-([0:N-1]*2*pi/N)));
    F_l = zeros(N, K);
    for i_k = 1:K
        ac_k = ac((i_k-1)*P+1:i_k*P);
        f_k = ac_k*exp(-1j*2*pi*r_radius/c*fc*cos(azimuth(i_m)-(pi*(2*i_k-1)/K)));
        F_l((i_k-1)*P+1:i_k*P,i_k) = f_k;
        t_k = r_radius/c*(1-cos(azimuth(i_m)-pi*(2*i_k-1)/K));
        t_k_list(i_k, i_m) = t_k;
    end
    F_PS = [F_PS, F_l];
end

end

