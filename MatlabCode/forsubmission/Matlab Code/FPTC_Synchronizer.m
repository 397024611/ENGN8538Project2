%% FPTC Synchronizer
% Description: Synchronize the timing and frequency offset. 
% The received TD OFDM sequence y and xp act as inputs and two outputs include
% the estimated timing n_est and the estimated frequency offset v_est.

function [n_est, v_est] = FPTC_Synchronizer(y, xp)
% Parameters
N = 2048; % FD OFDM symbol size
Ng = 512; %cyclic-prefix length
Dx = 12; % FD pilot distance
M = 4; % M samples in each segment
q = N/M; % divide the correlation into q segments

Nt = 100; 
Cor_FPTC = zeros(Ng + 2 * Nt + 1,q); % record data of segment correlation
Cor_FPTC_result = zeros(Ng + 2 * Nt + 1,1); % record data of the correlation with timings
% Entire y for correlation
for n = (N + Ng - Nt):(N + 2 * Ng + Nt) 
    for i = 1:q
        corr_seg = y((n + (i - 1) * M): (n + i * M - 1)) .*conj(xp(((i - 1) * M + 1): (i * M))); 
        Cor_FPTC(n - (N + Ng - Nt) + 1, i) = sum(corr_seg); % record data of segment correlation C(i)
    end
    Cor_FPTC_result(n - (N + Ng - Nt) + 1) = sum(abs(Cor_FPTC(n - (N + Ng - Nt) + 1, :)));
    % record the data of correlation with timings
end

[~, n_idx] = max(Cor_FPTC_result); % Obtain the estimated n
n_est = n_idx + N + Ng - Nt - 1; % Obtain the correct position for the end of 2nd CP


%% Fractional Frequency Synchronization (FFS) Using FPTC
% Method of Differential Correction
C = Cor_FPTC(n_idx, :); % get data of segment correlation for correct timings
Cor_diff = zeros(1, (q-1)); % record the data of differential correction
for j = 1:(q-1);
    Cor_diff(j) = C(j)*conj(C(j+1)); 
end
DC = mean(Cor_diff);

v_est = - N * angle(DC) / (2 * pi * M); % Obtain the estimated FFO v eventually

end
