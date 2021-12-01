%% Synthetic EEG dataset generation
% Mutiple source activity: Gaussian windowed sinusoidal activity

% Created by Andres Soler
% Reference biorXiv pre-print: "Automated methodology for optimal selection 
% of the minimum electrode subset for accurate EEG source localization 
% based on Genetic Algorithm", by Soler et al. (2021) 
% doi: https://doi.org/10.1101/2021.11.24.469917

% This code uses as forward model "The New York Head (ICBM-NY)" from "The 
% New York Head—A precise standardized volume conductor model for EEG 
% source localization and tES targeting" by Huang et al. (2015)
% doi: https://doi.org/10.1016/j.neuroimage.2015.12.019 available at
% https://www.parralab.org/nyhead/

% Variables in the output dataset:
% Fs: Sample frequency
% loc: Location of the source activity (voxel number)
% loc_xyz: Location of the source activity (3D space)
% t: Time
% Trials: Number of trials generated
% y: Simulated EEG


clear
close all
clc

%% EEG simulation Head Model with 231 Electrodes and 10K sources
%Forward head model, adapted from "the New York Head model" see ref above.
load('Head_Model_10K_231E.mat')

%% EEG Properties
%Sample Frequency (Hz)
Fs=200; 
%Simulation Interval (s)
inter=[0,3.5];
%Signal-to-noise SNR ratio (dB)
SNR = 0;
%Number of Trials
Trials = 150;
%Time 
t=inter(1)+(1/Fs):1/Fs:inter(2);

%% Source Activity Properties (Frequency, Time, Positions)
%Source 1
%Frequency s_1
f_s1=19;
%Center s_1
t_s1=0.5;
%Location s_1(Occipital Areas)
loc_s1 = [3727 8735 2734 7742 3461 8469 1205 6213 1180 6188 700 5708];
% Positions of loc_s1 verified that are in the 5k model

%Source 2
%Frequency s_2
f_s2=10;
%Center s_2
t_s2=1;
%Location s_2(Motor Cortex)
loc_s2 = [3422 8430 2284 7292 2271 7279 764 5772 1973 6981 266 5274];
% Positions of loc_s2 verified that are in the 5k model

%Source 3
%Frequency s_3
f_s3=7;
%Center s_3
t_s3=1.5;
%Location s_3(Frontal Cortex)
loc_s3 = [1784 6792 3511 8519 2153 7161 1600 6608 712 5720 3228 8236];
% Positions of loc_s3 verified that are in the 5k model

%Source 4
%Frequency s_4
f_s4=21;
%Center s_4
t_s4=2;
%Location s_4(Occipital Areas)
loc_s4 = loc_s1;

%Source 5 
%Frequency s_5
f_s5=12;
%Center s_5
t_s5=2.5;
%Location s_5(Motor Cortex)
loc_s5 = loc_s2;

%Source 6 (Frontal Cortex)
%Frequency s_6
f_s6=8;
%Center s_6
t_s6=3;
%Location s_6(Frontal Cortex)
loc_s6 = loc_s3;

%% Simulation Process
clc
disp('Starting EEG Simulation')
%EEG - initialization
y = zeros(size(M,1),size(t,2),Trials);
%Source activity - initialization
x = zeros(size(M,2),size(t,2),Trials);

%Random Amplitude
ub = 1; %upper bound
lb = 0.7; %lower bound
range = 10000;

for k=1:Trials
    
% Activity Source 1
temp1=exp(-0.5*(t-t_s1).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s1=temp1.*sin(2*pi*f_s1*t).*factor;
loc1 = loc_s1(ceil(rand*12));
x_s1=x__s1.*full(QG(:,loc1));
% Activity Source 2
temp2=exp(-0.5*(t-t_s2).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s2=temp2.*sin(2*pi*f_s2*t).*factor;
loc2 = loc_s2(ceil(rand*12));
x_s2=x__s2.*full(QG(:,loc2));
% Activity Source 3
temp3=exp(-0.5*(t-t_s3).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s3=temp3.*sin(2*pi*f_s3*t).*factor;
loc3 = loc_s3(ceil(rand*12));
x_s3=x__s3.*full(QG(:,loc3));
% Activity Source 4
temp4=exp(-0.5*(t-t_s4).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s4=temp4.*sin(2*pi*f_s4*t).*factor;
loc4 = loc_s4(ceil(rand*12));
x_s4=x__s4.*full(QG(:,loc4));
% Activity Source 5
temp5=exp(-0.5*(t-t_s5).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s5=temp5.*sin(2*pi*f_s5*t).*factor;
loc5 = loc_s5(ceil(rand*12));
x_s5=x__s5.*full(QG(:,loc5));
% Activity Source 6
temp6=exp(-0.5*(t-t_s6).^2/0.12^2);
factor =(lb+(ub-lb)*rand)/range;
x__s6=temp6.*sin(2*pi*f_s6*t).*factor;
loc6 = loc_s6(ceil(rand*12));
x_s6=x__s6.*full(QG(:,loc6));

% Save Source Positions
loc_xyz(:,:,k) = [vert(loc1,:);vert(loc2,:);vert(loc3,:);vert(loc4,:);vert(loc5,:);vert(loc6,:)];
loc(:,:,k) = [loc1; loc2; loc3; loc4; loc5; loc6];

% Source Activity
x(:,:,k) = x_s1+x_s2+x_s3+x_s4+x_s5+x_s6;

% EEG calculation
y(:,:,k)=awgn(M*x(:,:,k),SNR,'measured');
    
end
disp(['Data Simulated with SNR: ',num2str(SNR),' dB'])

%% Save Data
name = ['Dataset_SNR_0dB.mat'];
save(name,'y','loc_xyz','loc','t','Fs','Trials');
disp(['Saved File: ',name])