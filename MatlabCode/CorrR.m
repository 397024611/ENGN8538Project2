%Author : Yi Tang u5877586
%Date: 3 June 2017
%This function is used to do the correlation calculation which will be used in the MMSE.
function [R] = CorrR(x_Index,y_Index,Rho)
%"x_Index" is the index of one correlation input sequence.
%"y_Index" is the index of another correlation input sequence.
%"Rho" is the average power of each path.
N=1024;
N1=length(x_Index);
N2=length(y_Index);
Ld=length(Rho);
R=zeros(N1,N2);
for m=1:N1
    for n=1:N2
        R(m,n)=0;
        for l=0:1:Ld-1
            R(m,n)=R(m,n)+Rho(l+1)^2*exp(-2*j*pi*l*(x_Index(m)-y_Index(n))/N);
        end
    end
end
end

