%% FPTC Synchronizer
% Author: Hongli Shi u5738846
% Date: 10 May, 2017
% Description: The function is to synchronize the timing and frequency
% offset. The input is the received TD OFDM sequence y. There are two
% outputs. One is the estimated timing n_est (the end point of CP), 
% and the other is the estimated frequency offset v_est, including both  
% its integer and fractional parts.

function [n_est, v_est] = FPTC_Synchronizer(y, xp)
% y: the received TD signal from the channel
% xp: TD pilots
% n_est: the estimated timing offset
% v_est: the estimated frequency offset

% % Used for Test
% close all;
% clear all;
% clc

% Parameters
Dx = 6; % FD pilot distance
N = 1024; % FD OFDM symbol size, i.e., the FFt/IFFT size
Ng = 128; % the length of cyclic-prefix (CP)

% % Used for Test, Generate 3 Consecutive TD OFDM Symbols
% sym_num = 3; % the number of symbols generated
% [x_cp, xp, Xp, X_all] = OFDM_Symbol_Generator(sym_num);
% X = X_all(2,:);

% x_cp: 1 * sym_num * (N + Ng), the sym_num complete TD OFDM symbols with CP
% xp: 1 * N, TD pilots, same for all OFDM symbols
% Xp: 1 * N, FD pilots, same for all OFDM symbols
% X: sym_num * N, the 2nd FD OFDM symbol, used for test

%y = x_cp;

% % AWGN Channel
% SNR = 80;
% channel = 1; % Gaussian Channel
% [~,~,y] = ChannelSimulator(x_cp, SNR, channel, sym_num );
% %y = x_cp;
% 
% 
% % Frequency Offset
% v = 0;
% tt = 1:1:length(y); % time slot
% y = y .* exp(1i * 2 * pi * v * (tt-1) / N);

%% Coarse Timing Synchronization (CTS) Using FPTC
% Correlation Partition
M = 4; % each segment contains M samples
q = N / M; % partition the correlation to q segments
% Correlation = zeros((N + Ng),q); % record the segment correlation
% Correlation_result = zeros((N + Ng),1); % record the correlation under different timing
test_more = 100; % more correlation when starting point is before and after the 2nd CP
Correlation_FPTC = zeros(Ng + 2 * test_more + 1,q); % record the segment correlation
Correlation_FPTC_result = zeros(Ng + 2 * test_more + 1,1); % record the correlation under different timing

% for n = (N + Ng + 1):(2 * (N + Ng)) %
%     for i = 1:q
%         corr = y((n + (i - 1) * M): (n + i * M - 1)) .*...
%             conj(xp(((i - 1) * M + 1): (i * M))); % y_n_head[n]*conj(xp[n])
%         Correlation((n - N - Ng), i) = sum(corr); % record the segment correlation C(i)
%     end
%     Correlation_result(n - N - Ng) = sum(abs(Correlation((n - N - Ng), :)));
%     % record the correlation Corr under different timing 
% end

% Corr_result = []; % a vector to record the corr. under different timing
% n_est = 1001; % initialize our estimate on the CP starting location
% while (n_est+N+Ng-1 < 3200)
% corr = y(n_est: n_est+Ng-1) * y(n_est+N: n_est+Ng+N-1)';
% Corr_result = [Corr_result, corr];
% 
% n_est = n_est + 1; % increase the timing estimate by 1 and try again
% end
% [ignore, n_est_final] = max(abs(Corr_result)); % 
% figure(1)
% plot(1:1:length(Corr_result),abs(Corr_result));
% 
% % FFT for frequency domain signal vector
% y_est_0 = y(n_est_final+1000+Ng: n_est_final+1000+Ng+N-1);
% Y_0=fft(y_est_0, N);
% verify_0 = (sum(abs(Y_0-X))<0.00001); % Y should be equal to ifft input

for n = (N + Ng - test_more):(N + 2 * Ng + test_more) % correlate with the whole y
    for i = 1:q
        corr_seg = y((n + (i - 1) * M): (n + i * M - 1)) .*...
            conj(xp(((i - 1) * M + 1): (i * M))); % y_n_head[n]*conj(xp[n])
        Correlation_FPTC(n - (N + Ng - test_more) + 1, i) = sum(corr_seg); % record the segment correlation C(i)
    end
    Correlation_FPTC_result(n - (N + Ng - test_more) + 1) = sum(abs(Correlation_FPTC(n - (N + Ng - test_more) + 1, :)));
    % record the correlation Corr under different timing 
end

% figure(2)
% % Plot the Correlation under Different Timing
% plot(1:1:length(Correlation_FPTC_result), Correlation_FPTC_result); 

% Get Timing Offset
% n0 = N + 2 * Ng + 1; % the correct timing
[~, n_idx] = max(Correlation_FPTC_result); % find the estimated n that maximizes correlation
n_est = n_idx + N + Ng - test_more - 1; % find the true position for the end of 2nd CP
% n_diff = abs(n_est-n0);

% % Coarse Timing Synchronization (CTS)
% y_est_FPTC_CTS = y(n_est:(n_est+N-1)); % get the estimated TD y with correct timing
% Y_est_FPTC_CTS = fft(y_est_FPTC_CTS); % get the estimated FD Y with correct timing
% verify_1_FPCT_TS = (sum(abs(Y_est_FPTC_CTS - X))<0.00001); % verification: 1--Y_est = X; 0--Y_est ~= X.

%% Fractional Frequency Synchronization (FFS) Using FPTC
% Differential Correction
C = Correlation_FPTC(n_idx, :); % get the segment correlation for the correct timing
Differential_Correlation = zeros(1, (q-1)); % record the differential correction
for j = 1:(q-1);
    Differential_Correlation(j) = C(j)*conj(C(j+1)); % calculation according to Lecture 23
end
DC = mean(Differential_Correlation); % get the differential correction

% Estimated Frequency Offset Generation
% if(real(DC) < 0.00001) % ignore the real part which is too small
%     DC = 1i*imag(DC);
% end
% if(imag(DC) < 0.00001) % ignore the image part which is too small
%     DC = real(DC);
% end

v_est = - N * angle(DC) / (2 * pi * M); % get the estimated FFO v

end
