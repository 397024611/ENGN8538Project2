%Author : Yi Tang u5877586
%Date: 3 June 2017
%This function is the channel estimation function.
function [H_hat,Hd_hat] = ChannelEstimation(y,EstimationMethod,alpha,channel,SNR,n_hat)
%"y" is the received time domain signal sequence. 1*1152 complex double
%"EstimationMethod" is the selection of the estimation method. 1*1 double
%	1 : 1D-LS
%	2 : MMSE
%	3 : DTLS
%	4 : MMSE with SVD optimization
%"alpha" is the alpha coefficient which used in DTLS. 1*1 double
%"channel" is the selection of the channel type. 1*1 double
%	2 : Multi-path Rayleigh fading channel
%	3 : Multi-path Ricean fading channel
%	4 : Typical-urban 6-path channel
%"SNR" is the signal to noise ratio with unit 1. 1*1 double
%"n_hat" is the timing offset.

%Load from the data file to get the matrix parameter.
load('parameter.mat');
load('PilotPattern.mat');
%Parameters set up
N=1024;
Ng=128;
Dx=6;
beta=4/3;
pilot_position=[1:Dx:N,N];
Np=length(pilot_position);
data_position=setdiff(1:N,pilot_position);
Xp=PilotPattern;

%Calculate FD Y and estimation of pilots
Y=fft(y(n_hat+Ng+1:n_hat+Ng+1+N-1));
Hpls_hat=Y(pilot_position)./Xp;
H_hat=zeros(1,N);
H_hat(pilot_position)=Hpls_hat;

%Based on the channel type to choose Rhp and Rpp matrix from the data file.
%"P" is the number of paths with non-zero power coefficient in the channel. It will be used in the SVD optimization part.  
if channel==2
    Rhp=Rhp_ray;
    Rpp=Rpp_ray;
    P=12;
elseif channel==3
    Rhp=Rhp_Ricean;
    Rpp=Rpp_Ricean;
    P=13;
elseif channel==4
    Rhp=Rhp_Tu6;
    Rpp=Rpp_Tu6;
    P=6;
else
    H_hat=[];
    Hd_hat=[];
    warning('Channel Selection Error');
    return;
end

%Start estimation
if EstimationMethod==1
%1D_LS
    Hd_hat=interp1(pilot_position,Hpls_hat,data_position);
    H_hat(data_position)=Hd_hat;
elseif EstimationMethod==2
%MMSE
    H_hat=(Rhp*(Rpp+eye(Np)/(beta^2*SNR))^(-1)*Hpls_hat.').';
    Hd_hat=H_hat(data_position);
elseif EstimationMethod==3  
%DTLS
    h_hat=((FpHFp+alpha*eye(Np))^(-1)*FpH*Hpls_hat.').';
    H_hat=fft(h_hat(1:Ng),N);
    Hd_hat=H_hat(data_position);
elseif EstimationMethod==4
%MMSE with SVD optimization
    Rpphalf=(Rpp+eye(Np)/(beta^2*SNR))^(-0.5);
    RhpRpphalf=Rhp*Rpphalf;
    [U S V]=svd(RhpRpphalf);
	%Calculate part of the matrix and merge with zeros matrix to reduce the complexity.
	%Using S(1:P,1:P) to remove the extremely small and unnecessary values in the diagonal elements of the matrix S.
    H_hat=([U(:,1:P)*S(1:P,1:P),zeros(N,Np-P)]*V'*Rpphalf*Hpls_hat.').';
    Hd_hat=H_hat(data_position);
else
    H_hat=[];
    Hd_hat=[];
    warning('Estimation Method Selection Error');
end
end

