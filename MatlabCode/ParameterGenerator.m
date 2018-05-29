%Author : Yi Tang u5877586
%Date: 3 June 2017
%This script is a parameter generator. 
%It will generator the matrix which will be used in MMSE and DTLS then save the data in "parameter.mat".
%Once the file "parameter.mat" is generated, the task of this file is over. 
close all
clear all

%parameter set up
N=1024;
Ng=128;
Dx=6;
beta=4/3;
EXPTIMES=100;
pilot_position=[1:Dx:N,N];
Np=length(pilot_position);
x_cp=ones(1,N+Ng);

p=[0.248 0.129 0 0.31 0.425 0.49 0 0.0365 0.12 0 0 0 0.2 0 0 0 0 0.419 0 0 0 0 0 0 0.317 0 0 0 0 0.2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.185];
theta=[-2.57 -2.12 0 0.35 0.42 2.72 0 -1.44 1.13 0 0 0 -0.81 0 0 0 0 -1.55 0 0 0 0 0 0 -2.22 0 0 0 0 2.84 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.86];
% 'p','theta' for the fading channel in table 2.4 from A Study of DVB-T2 Standard with Physical Layer Transceiver Design and Implementation

%As the Rayleigh channel and the Ricean channel we use in this project is fixed, we can directly call the function "ChannelSimulator" to generate the channel coefficients.
[h_ray, ~, ~]=ChannelSimulator(x_cp, 1, 2, 1 );
Rho_ray=abs(h_ray);
%Call the "CorrR" function to calculate the Rhp and Rpp.
Rhp_ray=CorrR(1:N,pilot_position,Rho_ray);
Rpp_ray=CorrR(pilot_position,pilot_position,Rho_ray);

[h_Ricean, ~, ~]=ChannelSimulator(x_cp, 1, 3, 1 );
Rho_Ricean=abs(h_Ricean);
Rhp_Ricean=CorrR(1:N,pilot_position,Rho_Ricean);
Rpp_Ricean=CorrR(pilot_position,pilot_position,Rho_Ricean);

%For the Tu-6 channel, the coeffients are randomed.
%Here we simply call the "ChannelSimulator" function 100 times to calculate the average value as our expectation value of Rho.
for i=1:EXPTIMES
    [h_Tu6(i,:), ~, ~]=ChannelSimulator(x_cp, 1, 4, 1 );
    Rho(i,:)=abs(h_Tu6(i,:));
end
Rho_Tu6=sum(Rho)/100;
Rhp_Tu6=CorrR(1:N,pilot_position,Rho_Tu6);
Rpp_Tu6=CorrR(pilot_position,pilot_position,Rho_Tu6);


%Generator the Fourier operator which will be used in DTLS
for m=1:Np
    for n=1:Np
        Fp(m,n)=exp(-1j*pi*2*(pilot_position(m)-1)*(n-1)/N);
    end
end
FpH=Fp';
FpHFp=FpH*Fp;
save('parameter.mat','Rhp_ray','Rpp_ray','Rhp_Ricean','Rpp_Ricean','Rhp_Tu6','Rpp_Tu6','Fp','FpHFp','FpH');