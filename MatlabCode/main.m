%% simulation of CPcorrelation as 5.1

close all;
clear;
clc;


% Parameters
Dx = 12; % FD pilot distance
N = 2048; % FD OFDM symbol size, i.e., the FFt/IFFT size
Ng = 512; % the length of cyclic-prefix (CP)

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
plot(SNR, average_time_offset,  'r -*', 'LineWidth',2);
xlabel('SNR');
ylabel('Average timing offset  ');
set(gca,'FontSize',16);
grid on;
title('Average timing offset  ');


figure;
plot(SNR, num_ISI,  'r -*',  'LineWidth',2);
xlabel('SNR');
ylabel('Number of estimated timing');
grid on;
title('Estimated timing of incurring ISI');

figure;
semilogy(SNR, MSE,  'r -*',  'LineWidth',2);
xlabel('SNR');
ylabel('MSE');
grid on;
title('MSE of estimated  fractional frequency offset');

set(gca,'Fontsize', 20);
drawnow;
pause(0.1);



