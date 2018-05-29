%Author : Yi Tang u5877586
%Date: 3 June 2017
%This script is used to generate the pilot pattern.
close all
clear all

%In our project, we make the random seed equals to 1 to make sure that every time we can generate same pilot sequence.
rand('seed',1);
beta=4/3;
PilotPattern = beta*(randi([0, 1], 1, 172)*2-1);
save ('PilotPattern.mat','PilotPattern');