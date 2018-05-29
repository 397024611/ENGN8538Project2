%% Description
% % y1: the received TD signal for 1symbol
% % y2: the received TD signal for 2symbols
% % Xp1: FD pilots for 1symbol
% % Xp2: FD pilots for 2symbols
% % TO1: the timing offset for 1symbol
% % TO2: the timing offset for 2symbols
% % vi_est_1sym: the estimated integer frequency offset for 1symbol
% % vi_est_2sym: the estimated integer frequency offset for 2symbols
%% function
function [vi_est_1sym, vi_est_2sym] = Integer_Frequency_Synchronizer(y1, y2, Xp1, Xp2, TO1, TO2)

Dx = 12;
N = 2048;
Ng = 512;

n0 = N + 2 * Ng + 1; % the correct timing

y_vi_2_1sym = y1(n0-TO1:(n0-TO1 + N - 1)); % the 2nd TD OFDM

y_vi_2_2sym = y2(n0-TO2:(n0-TO2 + N - 1)); % the 2nd TD OFDM
y_vi_3_2sym = y2((n0-TO2 + N + Ng):(n0-TO2 + 2 * N + Ng - 1)); % the 3rd TD OFDM

Y_2_1sym = fft(y_vi_2_1sym,N); % the 2nd FD OFDM for IFS using 1symbol
Y_2_2sym = fft(y_vi_2_2sym,N); % the 2nd FD OFDM for IFS using 2symbols
Y_3_2sym = fft(y_vi_3_2sym,N); % the 3rd FD OFDM for IFS using 2symbols

%% IFS Using 1Symbol
f_index = 1:1:N; % frequency index
Corr_1sym = zeros(N,1); %empty cell use to record the correlation result using 1symbol

for vi = 1:N
    Corr_1sym(vi) = sum(circshift(Y_2_1sym,(vi-1),2) .* conj(Xp1(f_index)));  
end
[~, vi_idx_1] = max(abs(Corr_1sym)); % find the estimated vi wiht max |Corr|
vi_est_1sym = vi_idx_1 - 1;
%% IFS Using 2Symbols
Z = Y_2_2sym .* conj(Y_3_2sym); % calculate Z
Corr_2 = zeros(N,1); %empty cell use to record the correlation result using 2symbols
for vi = 1:N
    Corr_2(vi) = sum(circshift(Z,(vi-1),2) .* conj((Xp2(f_index)).^2));
end
[~, vi_idx_2] = max(abs(Corr_2)); % same as 1symbol
vi_est_2sym = vi_idx_2 - 1;
end

