%% Integer Frequency Synchronizer (IFS)
% Author: Hongli Shi u5738846
% Date: 31 May, 2017
% Description: The function is to synchronize the integer frequency.

function [vi_est_1, vi_est_2] = Integer_Frequency_Synchronizer(y, y2, Xp, Xp2, TO_1, TO_2)
% y: the received TD signal from the channel (without continual pilots)
% y2: the received TD signal from the channel (with continual pilots)
% Xp: FD pilots (without continual pilots)
% Xp2: FD pilots (with continual pilots)
% TO_1: the timing offset for IFS using 1 symbol
% TO_2: the timing offset for IFS using 2 symbols
% vi_est_1: the estimated integer frequency offset for 1 symbol
% vi_est_2: the estimated integer frequency offset for 2 symbols

% Parameters
Dx = 6; % FD pilot distance
N = 1024; % FD OFDM symbol size, i.e., the FFt/IFFT size
Ng = 128; % the length of cyclic-prefix (CP)

% Correct Timing
n0 = N + 2 * Ng + 1; % the correct timing

% % Correct FFO
% % v_f = v - fix(v); % the correct fractional frequency offset
% TO = 50;

% Record y with Corresponding Timing Offset
y_vi_2_1 = y(n0-TO_1:(n0-TO_1 + N - 1)); % the 2nd TD OFDM

y_vi_2 = y2(n0-TO_2:(n0-TO_2 + N - 1)); % the 2nd TD OFDM
y_vi_3 = y2((n0-TO_2 + N + Ng):(n0-TO_2 + 2 * N + Ng - 1)); % the 3rd TD OFDM

% % Remove Fractional Frequency Offset
% t = 1:1:N; % time index
% y_vi_2 = y_n0_2 .* exp(-1i * 2 * pi * v_f * (t-1) / N); % the 2nd TD OFDM
% y_vi_3 = y_n0_3 .* exp(-1i * 2 * pi * v_f * (t-1) / N); % the 3rd TD OFDM

% Get FD Y with Only IFO through FFT
Y_2_1 = fft(y_vi_2_1,N); % the 2nd FD OFDM for IFS using 1 symbol
Y_2 = fft(y_vi_2,N); % the 2nd FD OFDM for IFS using 2 symbols
Y_3 = fft(y_vi_3,N); % the 3rd FD OFDM for IFS using 2 symbols

% % Verify Y_2 = X
% verify_2_IFS_FFS_nest = (sum(abs(Y_2-X))<0.00001);

%% Integer Frequency Synchronization (IFS) Using 1 Symbol
k = 1:1:N; % frequency index
Corr_1 = zeros(N,1); % record the correlation result using 1 symbol

% Calculate the Correlation
for vi = 1:N
    
%     if(vi==1)
%         Corr_1(vi) = sum(Y_2(k) .* conj(Xp(k)));
%     else
%         Corr_1(vi) = sum([Y_2((end - vi + 2):(end)),Y_2(1:(end - vi + 1))] .* conj(Xp(k)));
%     end
    
    Corr_1(vi) = sum(circshift(Y_2_1,(vi-1),2) .* conj(Xp(k)));
    
end

% figure(3)
% % Plot the Absulote Values of Correlation under Different IFO
% plot(1:1:length(Corr_1), abs(Corr_1));

% Get IFO
[~, vi_idx_1] = max(abs(Corr_1)); % find the estimated vi that maximizes |Corr|
vi_est_1 = vi_idx_1 - 1;

% if(vi_est_1 == 0)
%     Y_est_IFS_1 = Y_2;
% else
%     Y_est_IFS_1 = [Y_2((end - vi_est_1 + 1):(end)),Y_2(1:(end - vi_est_1))];
% end

% verify_3_IFS_1 = (sum(abs(Y_est_IFS_1 - X)) < 0.00001);

%% Integer Frequency Synchronization (IFS) Using 2 Symbols
Z = Y_2 .* conj(Y_3); % remove timing offset estimation error
Corr_2 = zeros(N,1); % record the correlation result using 2 symbols
% Note: 

% Calculate the Correlation
for vi = 1:N
    
%     if(vi==1)
%         Corr_2(vi) = sum(Z(k) .* conj(Xp(k)));
%     else
%         Corr_2(vi) = sum([Z((end - vi + 2):(end)),Z(1:(end - vi + 1))] .* conj(Xp(k)));
%     end
    Corr_2(vi) = sum(circshift(Z,(vi-1),2) .* conj(Xp2(k)));
end

% figure(4)
% % Plot the Absulote Values of Correlation under Different IFO
% plot(1:1:length(Corr_2), abs(Corr_2));

% Get IFO
[~, vi_idx_2] = max(abs(Corr_2)); % find the estimated vi that maximizes |Corr|

vi_est_2 = vi_idx_2 - 1;

% if(vi_est_2 == 0)
%     Y_est_IFS_2 = Y_2;
% else
%     Y_est_IFS_2 = [Y_2((end - vi_est_2 + 1):(end)),Y_2(1:(end - vi_est_2))];
% end

% verify_4_IFS_2 = (sum(abs(Y_est_IFS_2 - X)) < 0.00001);
end

