%Author : Yi Tang u5877586
%Date: 3 June 2017
%This function is used to generate the OFDM symbol with continuous pilots, 
%which will be used in IFS using two OFDM Symbols.
function [x_cp, Xp] = SymGen_ContinuousPilot(M)
%"M" is the number of the symbols we want to generate. 1*1 double.
%"x_cp" is the time domain sequence of the signal. M*1152 complex double.

N=1024;
Ng=128;%the length of the cyclic-prefix
Dx = 6; % pilot distance
beta=4/3;
scattered_pilot_carrier = [1: Dx: N,N];
temp_carrier = setdiff(1: N, scattered_pilot_carrier);
randomIndex=randperm(length(temp_carrier));
continual_pilot_carrier = temp_carrier(randomIndex(1:20));
pilot_carrier=[scattered_pilot_carrier continual_pilot_carrier];
Np=length(pilot_carrier);
Xp=zeros(1,N);
Xp(pilot_carrier)=beta;
data_carrier = setdiff(1: N, pilot_carrier);
%% Generate the data from QPSK
data_I=(randi([0, 1], M, length(data_carrier))*2-1)*(sqrt(2)/2);
data_Q=(randi([0, 1], M, length(data_carrier))*2-1)*(sqrt(2)/2);
%% generate the frequency domain signal vector
X = zeros(M, N);
X(:,data_carrier) = data_I + 1i*data_Q;
%% generte the pilot pattern
X(:,pilot_carrier) = repmat(beta*(1:Np),M,1); % boosted pilot
for i=1:M
%% IFFT for time domain signal vector
x(i,:)= ifft(X(i,:), N);
% generate cyclic-prefix attached time domain signal vector
xcp(i,:) = [x(i,end-Ng+1: end), x(i,:)];
end
xcp=[x(:,N-Ng+1:N) x];
x_cp=reshape(xcp.',1,numel(xcp));
end
