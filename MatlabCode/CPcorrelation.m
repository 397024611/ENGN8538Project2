%Author : Chuan Qin, u5832845
function [delta_n,v_f]=CPcorrelation(v,SNR,M)
% 'v'    frequency offset
% 'M'    the number of OFDM symbols 
% 'SNR'  is in dB
% 'v_f'  fractional frequency offset
% 'delta_n'   the difference between the tested and theoretical value
N=2048;
Ng=512;
Dx = 12;          % piolt distance
M=3;
   x_cp=Symbol_generator(M);                  % though the symbol generator
   x_cp_s = reshape(x_cp.',1,[]);     
   y=awgn(conv(1,x_cp_s), SNR, 'measured');   % though the AWGN channel
   y_1=reshape(y,[],M);
   x_cp_index=repmat(1:N+Ng,M,1);             % make 3 x_cp

y_2=(y_1.').*exp(1i*2*pi*v*x_cp_index/N);     %timing synchronization, IFO+FFO+NOISE
y=[y_2(1,end-Ng-1:end),y_2(2,:),y_2(3,1:Ng)]; % 3 times of Ng

Corr_result = [];                             % a vector to record the corr. under different timing
n_est = 1;                                    % initialize our estimate on the CP starting location

while (n_est+N+Ng-1<length(y))
corr=y(n_est:n_est+Ng-1)*y(n_est+N: n_est+Ng+N-1)';
Corr_result=[Corr_result, abs(corr)];

n_est = n_est + 1;                            % increase the timing estimate by 1 and try again
end
% figure;
% stem(N+Ng+1:N+Ng+length(Corr_result), Corr_result);
% xlim([0, 2*(N+Ng) + 1]); ylim([0, 0.15]);
[~, n_est_final] = max(Corr_result);  

corr_final = y(n_est_final: n_est_final+Ng-1)* y(n_est_final+N: n_est_final+Ng+N-1)';
C=angle(corr_final)/(-2*pi);
if C<0
    C=1-abs(C);
else
    C=C;
end
v_f=C;
delta_n=abs(n_est_final-Ng);
