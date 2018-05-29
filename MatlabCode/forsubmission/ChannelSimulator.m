% 'x_cp' TD OFDM symbol with CP
% 'M'    the number of OFDM symbols 
% 'SNR'  is in dB
% 'h'    TD channel CIR
% 'Ld'   the length of h
% 'y'    the received TD OFDN symbol


function [h, Ld, y]=ChannelSimulator(x_cp, SNR, channel, M )

h_ray=zeros(1,50);
h_ray([1,2,4,5,6,8,9,13,18,25,30,50]) = [0.248*exp(2.57i),0.129*exp(2.12i),...
    0.13*exp(-0.35i),0.425*exp(-0.42i),0.49*exp(-2.72i),0.0365*exp(1.44i),...
    0.12*exp(-1.13i),0.2*exp(0.81i),0.419*exp(1.55i),0.317*exp(-2.22i),...
    0.2*exp(-2.84i),0.185*exp(-2.86i)];   % from table 2.4 of standard  
h_ray=h_ray./sqrt(sum(abs(h_ray).^2));    % normalize the total power to 1

if channel==1                         % AWGN channel
    h=1;
    Ld=1;
    y=awgn(conv(h, x_cp), SNR, 'measured');
    
elseif channel==2                     % Multi-path Rayleigh fading channel
    Ld=50;
    h=h_ray;
    y=awgn(conv(h, x_cp), SNR, 'measured'); %convolute the TD OFDM symbol with h and add noise
    
elseif channel==3                     % Multi-path Ricean fading channel
    h=[0 h_ray];
    h(1)=sqrt(10*sum(abs(h_ray).^2));   % the first tap of h
    h=h./sqrt(sum(abs(h).^2));
    Ld=51;
    y=awgn(conv(h, x_cp), SNR, 'measured');  %convolute the TD OFDM symbol with h and add noise
else
    error('Error channel type');
end
end


