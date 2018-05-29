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
