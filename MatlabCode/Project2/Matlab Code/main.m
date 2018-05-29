% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function for project two, which runs the simulation of OFDM
% symbols by using three synchronization techniques
% author: Peng Chen
%         Tuo Zhao
%         Jiawei Li
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;
% Parameters used for this project (group K)
Dx = 12; % FD pilot distance
N = 2048; % FD OFDM symbol size, i.e., the FFt/IFFT s ize
Ng = 512; % the length of cyclic-prefix (CP)

%% simulation of CPcorrelation as 5.1
n_act=N+Ng+1;
v=4.9;
vf_act=v-floor(v);
sim_number=100;
SNR=0:3:30;

average_time_offset=zeros(1, length(SNR));
MSE=zeros(1, length(SNR));
delta_n=zeros(1, sim_number);
v_f=zeros(1, sim_number);
num_ISI=zeros(1, length(SNR));
for j=1:length(SNR)
    for i=1:sim_number
        [delta_n(i), v_f(i)]=CPcorrelation(v, SNR(j),3);
    end
    ISI=(delta_n+N+Ng)-n_act;
    num_ISI(j)=length(find(ISI>0));
    average_time_offset(j)=mean(abs(ISI));
    MSE(j)=sum((v_f-vf_act).^2)/sim_number;
end

figure;
plot(SNR, average_time_offset,  '-*', 'LineWidth',2);
xlabel('SNR');
ylabel('Average timing offset  ');
title('Average timing offset  ');


figure;
plot(SNR, num_ISI,  '-*',  'LineWidth',2);
xlabel('SNR');
ylabel('Number of estimated timing');
title('Estimated timing of incurring ISI');

figure;
semilogy(SNR, MSE,  '-*',  'LineWidth',2);
xlabel('SNR');
ylabel('MSE');
title('MSE of estimated  fractional frequency offset');
drawnow;
pause(0.1);

%% 5.2
clear;
clc;

% Parameters
N = 2048; % FD OFDM symbol size
Ng = 512; %cyclic-prefix length
Dx = 12; % FD pilot distance
no_sym = 3; % Generate 3 TD OFDM Symbols
[x_cp, xp, Xp] = OFDM_Symbol_Generator(no_sym);

% AWGN Channel parameters
SNR = 0:3:30;
channel = 1; 

% Frequency Offset parameters
v = 2.5;
ts = 1:1:length(x_cp);

% Correct Timing
n0 = N + 2 * Ng + 1; % The correct timing
sim_times = 20; % Run 10 times 


avg_TO = zeros(1, length(SNR)); % Get average timing offset
MSE_FO = zeros(1, length(SNR)); % Get the MSE of estimated frequency offset

for j = 1:length(SNR)
    % Obtain received signal y from AWGN channel with different SNRs
    [~,~,y] = ChannelSimulator(x_cp, SNR(j), channel, no_sym );
    
    % Determine a frequency offset for received signal y
    y = y .* exp(1i * 2 * pi * v * (ts-1) / N);
    
    n_diff = zeros(1, sim_times); % Time offset
    v_diff = zeros(1, sim_times); 
    
    % Get the square of estimation error of FO
    for i = 1:sim_times
        [n_est, v_est] = FPTC_Synchronizer(y, xp); % FPTC
        % Get Timing Offset
        n_diff(i) = abs(n_est-n0); % The estimation error of time
        v_diff(i) = (abs(v_est-v))^2; % The square estimation error of frequency
    end
    avg_TO(j) = mean(n_diff, 2); % Average timing offset
    MSE_FO(j) = mean(v_diff, 2); % The MSE of estimated frequency offset
end

figure;
plot(SNR, avg_TO,  '-*', 'LineWidth',2);
xlabel('SNR/dB');
ylabel('average timing offset');
title('Average Timing Offset v.s. SNR');
set(gca,'FontSize',14);
grid on;


figure;
plot(SNR, MSE_FO,  '-*',  'LineWidth',2);
xlabel('SNR/dB');
ylabel('MSE of estimated frequency offset');
title('MSE of Estimated Frequency Offset v^(^) v.s. SNR');
set(gca,'FontSize',14);
grid on;


set(gca,'Fontsize', 20);
drawnow;
pause(0.1);


%% 5.3
clear;
clc;
Dx = 12; % FD pilot distance
N = 2048; % FD OFDM symbol size, i.e., the FFt/IFFT size
Ng = 512; % the length of cyclic-prefix (CP)
% Setting Timing Offset
TO1 = 0; % timing offset for IFS using 1 symbol, no TO for 1 symbol
TO2 = 5; % timing offset for IFS using 2 symbols
% Generate 3 Consecutive TD OFDM Symbols
sym_num = 3; % the number of symbols generated
[x_cp1,~, Xp1] = OFDM_Symbol_Generator(sym_num);
[x_cp2,~, Xp2] = OFDM_Symbol_Generator(sym_num);
% AWGN Channel
SNR = 0:2:30; % SNR setting in db
channel = 1; % Gaussian Channel
% Frequency Offset
v = 10;
tts = 1:1:length(x_cp1); % time slot
% Correct Timing
n0 = N + 2 * Ng + 1; % the correct timing
sim_times = 10; % The number times for simulation

MSE_vi_1 = zeros(1, length(SNR)); % record the MSE of estimated integer frequency offset for 1 symbol
MSE_vi_2 = zeros(1, length(SNR)); % record the MSE of estimated integer frequency offset for 2 symbols

for j = 1:length(SNR)
    [~,~,y1] = ChannelSimulator(x_cp1, SNR(j), channel, sym_num );
    [~,~,y2] = ChannelSimulator(x_cp2, SNR(j), channel, sym_num );
    y1 = y1 .* exp(1i * 2 * pi * (-v) * (tts-1) / N);
    y2 = y2 .* exp(1i * 2 * pi * (-v) * (tts-1) / N);
    
    vi_diff_1sym = zeros(1, sim_times); 
    vi_diff_2sym = zeros(1, sim_times); 
    for i = 1:sim_times
        [vi_est_1sym, vi_est_2sym] = Integer_Frequency_Synchronizer(y1, y2, Xp1, Xp2, TO1, TO2);
        vi_diff_1sym(i) = (abs(vi_est_1sym-v))^2; %square estimation error
        vi_diff_2sym(i) = (abs(vi_est_2sym-v))^2; %square estimation error 
    end
    MSE_vi_1(j) = mean(vi_diff_1sym, 2); % record MSE for 1 symbol
    MSE_vi_2(j) = mean(vi_diff_2sym, 2); % record MSE for 2 symbols
end
%
    figure;
    plot(SNR, MSE_vi_1,  'b -o', 'LineWidth',1.8);
    title('MSE of Estimated Integer Frequency Offset VS SNR(1Symbol)');
    xlabel('SNR(dB)');
    ylabel('MSE of estimated integer frequency offset');
    grid on;
%
    figure;
    plot(SNR, MSE_vi_2,  'm -o',  'LineWidth',1.8);
    title('MSE of Estimated Integer Frequency Offset VS SNR (2Symbols)');
    xlabel('SNR(dB)');
    ylabel('MSE of estimated integer frequency offset');
    grid on;

