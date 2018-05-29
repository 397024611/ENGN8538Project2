%Author : Yi Tang u5877586
%		  Chuan Qin u5832845 
%Date: 3 June 2017
%This function is used to generate the OFDM symbol.
function x_cp=Symbol_generator(M)
%"M" is the number of the symbols we want to generate. 1*1 double.
%"x_cp" is the time domain sequence of the signal. M*1152 complex double.
load('PilotPattern.mat');
N=2048 ;
Ng=512;%the length of the cyclic-prefix
Dx = 12; % pilot distance

%% generate the location of pilot carriers
pilot_carrier = [1: Dx: N,N]; % edge pilots are included
%% generate the location of data carriers
data_carrier = setdiff(1: N, pilot_carrier);
Np=length(pilot_carrier);

%% Generate the data from QPSK
data_I=(randi([0, 1], M, length(data_carrier))*2-1)*(sqrt(2)/2);
data_Q=(randi([0, 1], M, length(data_carrier))*2-1)*(sqrt(2)/2);
%% generate the frequency domain signal vector
X = zeros(M, N);
X(:,data_carrier) = data_I + 1i*data_Q;
%% generte the pilot pattern
X(:,pilot_carrier) = repmat(PilotPattern,M,1); % boosted pilot
for i=1:M
%% IFFT for time domain signal vector
x(i,:)= ifft(X(i,:), N);
% generate cyclic-prefix attached time domain signal vector
x_cp(i,:) = [x(i,end-Ng+1: end), x(i,:)];
end
x_cp=[x(:,N-Ng+1:N) x];
end