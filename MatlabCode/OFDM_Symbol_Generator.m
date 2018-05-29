%% OFDM Symbol Generator
% Author: Hongli Shi u5738846
% Date: 10 May, 2017
% Description: The function is to generate the complete sym_num TD OFDM 
% symbols xcp with cyclic prefix. An FD OFDM symbol is first generated, 
% which contains two types of carriers: pilot (known by the receiver, same 
% pilots for all OFDM symbols) and data. FD OFDM symbol X is directly 
% IFFT-ed to obtain the TD symbol x. The last Ng samples of x is copied and 
% attached to the front of x as CP. The output of the function is sym_num 
% complete TD OFDM symbols x_cp, whose length is sym_num * (N + Ng).

function [x_cp, xp, Xp] = OFDM_Symbol_Generator(sym_num)

% sym_num: the number of symbols generated
% x_cp: 1 * sym_num * (N + Ng), the sym_num complete TD OFDM symbols with CP
% xp: 1 * N, TD pilots, same for all OFDM symbols
% Xp: 1 * N, FD pilots, same for all OFDM symbols

% sym_num = 3; % used for test

% Parameters
Dx = 6; % FD pilot distance
N = 1024; % FD OFDM symbol size, i.e., the FFt/IFFT size
Ng = 128; % the length of cyclic-prefix (CP)

% FD OFDM Symbol Generation
X = zeros(1, N); % record one FD OFDM symbol%
% X_all = zeros(sym_num, N); % record all the FD OFDM symbols
QPSK_table = [ (sqrt(2)/2) + 1i * (sqrt(2)/2), (sqrt(2)/2) - 1i * (sqrt(2)/2),...
    - (sqrt(2)/2) + 1i * (sqrt(2)/2), - (sqrt(2)/2) - 1i * (sqrt(2)/2)]; % DVB-T2
% QPSK table stores the complex values for the four cells.

% Positions for Pilots and Data
pilot_position = [1:Dx:N,N]; % get the positions for pilots
data_position = setdiff(1:N,pilot_position); % get the positions for data

% Read pilot pattern from PilotPattern.mat.
load('PilotPattern.mat');
% the amplitude of all pilots is 4/3, each pilot has different signs.

% Pilots Assignment
X(pilot_position) = PilotPattern; % assign the edge pilots and scattered pilots
Xp = X; % get the FD pilots with all data = 0

% sym_num Complete TD OFDM Symbols Generation
x_cp = zeros(1, sym_num * (N + Ng)); % record the sym_num TD OFDM symbols
% x_all = zeros(sym_num, N); % record the sym_num TD OFDM symbols without CP
for i = 1:sym_num
% Data Generation
% time = clock; % 
% rand('seed', time(6)); % use the value of second as the seed 
data_idx = randi([1,4],1, length(data_position)); % generate random numbers 1, 2, 3, 4.
data = QPSK_table(data_idx); % get complex values from the QPSK table.

% Data Assignment
X(data_position) = data; % assign the data
% X_all(i,:) = X; % record the i-th FD OFDM symbol to the i-th row

% TD OFDM Symbol Generation
x = ifft(X, N); % get TD OFDM symbol
% x_all(i,:) = x;
x_cp((i - 1) * (N + Ng) + 1: i * (N + Ng)) = [x((N-Ng+1):end), x]; 
% get the complete TD OFDM symbol by putting Cyclic Prefix ahead
end

%X2 = X_all(2,:); % record the 2nd FD OFDM symbol

% TD Pilots xp Generation
xp = zeros(1, N); % record TD pilots
for n = 1:N % time index
    xp(n) = (1/(N)) * sum((X(pilot_position) .* exp(1i*2*pi*(n-1)*(pilot_position-1)/N)),2); % TD pilots
end

% Verify the Relationship x = xp + xd
% TD Data Generation
% tt = 1:1:N; % frequency index
% xd = zeros(sym_num, N); % record TD data
% xt = zeros(sym_num, N); % record TD OFDM symbol, i.e., ifft of FD OFDM symbol
% verification = zeros(sym_num,1);
% for j = 1:sym_num
%     for n = 1:N % time index
%         xd(j,n) = (1/(N)) * sum((X_all(j,data_position) .* exp(1i*2*pi*(n-1)*(data_position-1)/N)),2); % TD data
%         xt(j,n) = (1/(N)) * sum((X_all(j,tt) .* exp(1i*2*pi*(n-1)*(tt-1)/N)),2); % TD OFDM symbol
%     end
%     verification(j) = (sum(abs(xp+xd(j,:)-x_all(j,:)),2)<0.00001); % verification: 1--x = xp + xd; 0: x ~= xp + xd
% end
% 
end
