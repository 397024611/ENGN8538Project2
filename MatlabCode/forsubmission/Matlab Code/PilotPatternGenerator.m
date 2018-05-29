
close all
clear all

rand('seed',1);
% Set the amplitude to 4/3
beta=4/3;
PilotPattern = beta*(randi([0, 1], 1, 172)*2-1);
save ('PilotPattern.mat','PilotPattern');