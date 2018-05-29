% % sym_num: How many OFDM continuous samples will used to generate
% % x_cp: 1 * sym_num * (N + Ng), the sym_num complete TD OFDM symbols with CP
% % xp: 1 * N, TD pilots
% % Xp: 1 * N, FD pilots
% 
% % sym_num = 3 in this project
% function [x_cp, xp, Xp] = OFDM_Symbol_Generator(sym_num)
% 
% 
% % %% Parameters
% Dx = 12;  % FD pilot distance
% N = 2048; % FD OFDM symbol size, i.e., the FFt/IFFT size
% Ng = 512; % the length of cyclic-prefix (CP)
% 
% % FD OFDM Symbol Generation
% X = zeros(1, N); % record one FD OFDM symbol%
% % X_all = zeros(sym_num, N); % record all the FD OFDM symbols
% QPSK_table = [ (sqrt(2)/2) + 1i * (sqrt(2)/2), (sqrt(2)/2) - 1i * (sqrt(2)/2),...
%     - (sqrt(2)/2) + 1i * (sqrt(2)/2), - (sqrt(2)/2) - 1i * (sqrt(2)/2)]; % DVB-T2
% % QPSK table stores the complex values for the four cells.
% 
% % Positions for Pilots and Data
% pilot_position = [1:Dx:N,N]; % get the positions for pilots
% data_position = setdiff(1:N,pilot_position); % get the positions for data
% 
% % Read pilot pattern from PilotPattern.mat.
% load('PilotPattern.mat');
% % the amplitude of all pilots is 4/3, each pilot has different signs.
% 
% % Pilots Assignment
% X(pilot_position) = PilotPattern; % assign the edge pilots and scattered pilots
% Xp = X; % get the FD pilots with all data = 0
% 
% % sym_num Complete TD OFDM Symbols Generation
% x_cp = zeros(1, sym_num * (N + Ng)); % record the sym_num TD OFDM symbols
% % x_all = zeros(sym_num, N); % record the sym_num TD OFDM symbols without CP
% for i = 1:sym_num
% % Data Generation
% % time = clock; % 
% % rand('seed', time(6)); % use the value of second as the seed 
% data_idx = randi([1,4],1, length(data_position)); % generate random numbers 1, 2, 3, 4.
% data = QPSK_table(data_idx); % get complex values from the QPSK table.
% 
% % Data Assignment
% X(data_position) = data; % assign the data
% % X_all(i,:) = X; % record the i-th FD OFDM symbol to the i-th row
% 
% % TD OFDM Symbol Generation
% x = ifft(X, N); % get TD OFDM symbol
% % x_all(i,:) = x;
% x_cp((i - 1) * (N + Ng) + 1: i * (N + Ng)) = [x((N-Ng+1):end), x]; 
% % get the complete TD OFDM symbol by putting Cyclic Prefix ahead
% end
% 
% %X2 = X_all(2,:); % record the 2nd FD OFDM symbol
% 
% % TD Pilots xp Generation
% xp = zeros(1, N); % record TD pilots
% 
% for n = 1:N % time index
%     xp(n) = (1/(N)) * sum((X(pilot_position) .* exp(1i*2*pi*(n-1)*(pilot_position-1)/N)),2); % TD pilots
% end
% 
% end
%%
% symbol_num: the number of symbols generated
% x_cp: the sym_num complete TD OFDM symbols with CP
% xp: TD pilots, same for all OFDM symbols
% Xp: FD pilots, same for all OFDM symbols
%% Function OFDM Symbol Generator
function [x_cp, xp_TD, xp_FD] = OFDM_Symbol_Generator(symbol_num)

    Dx = 12; 
    N = 2048;
    Ng = 512; 

    X = zeros(1, N);
    QPSK_table = [ (sqrt(2)/2) + 1i * (sqrt(2)/2),... % QPSK table stores the complex values for the four cells.
        (sqrt(2)/2) - 1i * (sqrt(2)/2),...
        - (sqrt(2)/2) + 1i * (sqrt(2)/2), - (sqrt(2)/2)...
        - 1i * (sqrt(2)/2)]; % DVB-T2

    pilot_position = [1:Dx:N,N]; % get the positions for pilots

    data_position = setdiff(1:N,pilot_position); % get the positions for data

    load('PilotPattern.mat');% amplitude = 4/3, pilot may has different signs.

    X(pilot_position) = PilotPattern; % assign the edge pilots, scattered pilots

    xp_FD = X; % get the FD pilots with all data = 0

    x_cp = zeros(1, symbol_num * (N + Ng));

    for i = 1:symbol_num
    data_idx = randi([1,4],1, length(data_position)); % generate random integer 1 to 4
    data = QPSK_table(data_idx); % complex values from QPSK table.
    X(data_position) = data; % assign the data
    x = ifft(X, N); %TD OFDM symbol
    x_cp((i - 1) * (N + Ng) + 1: i * (N + Ng)) = [x((N-Ng+1):end), x]; % putting Cyclic Prefix ahead
    end

    xp_TD = zeros(1, N);

    for n = 1:N
        xp_TD(n) = (1/(N)) * sum((X(pilot_position) .* exp(1i*2*pi*(n-1)*(pilot_position-1)/N)),2); % TD pilots
    end
end

