% %% input:  v    frequency offset
% %%         M    the number of OFDM symbols 
% %%        SNR  is in dB
% %% output: et: estimated timing
% %%         ef: estimated frequency offeset 

function [et,ef]=CPcorrectelation(v,SNR,M)
N=2048;           % FD OFDM symbol size
Ng=512;           % the length of cyclic-prefix (CP)
Dx = 12;          % FD pilot distance
M=3;

% Generate OFDM symbol and pass through AMGN channel
[x_cp, ~, ~]=OFDM_Symbol_Generator(M);   
x_cp_s = reshape(x_cp.',1,[]);   
[~, ~, y] = ChannelSimulator(x_cp_s, SNR, 1, M );

y_1=reshape(y,[],M);                           %size : 2560 * 3
x_cp_index=repmat(1:N+Ng,M,1);                 % indicate the index of each x_cp symbol
y_ts=(y_1.').*exp(1i*2*pi*v*x_cp_index/N);     % from equation in lecture 23,p4, timing synchronization
y=[y_ts(1,end-(Ng-1):end),y_ts(2,:),y_ts(3,1:Ng)]; % recombine the symbol to cyclix prefix

result_array = [];                            
Estmated_n = 1;                                    

while (Estmated_n+N+Ng-1<length(y))
correct=y(Estmated_n:Estmated_n+Ng-1)*y(Estmated_n+N: Estmated_n+Ng+N-1)';
result_array=[result_array, abs(correct)];

Estmated_n = Estmated_n + 1;                            % increase the timing estimate by 1 and try again
end

[~, Estmated_n_final] = max(result_array);  

correct_final = y(Estmated_n_final: Estmated_n_final+(Ng-1)) * y(Estmated_n_final+N: Estmated_n_final+N+(Ng-1))';
C=angle(correct_final)/(-2*pi);
if C<0
    C=1-abs(C);
else
    C=C;
end
ef=C;
et=abs(Estmated_n_final-Ng);
