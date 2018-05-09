clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Separating into gait cycles'));
addpath(genpath('Healthy Recordings\Subject6\NO_FLOAT'));
addpath(fullfile('Separating into gait cycles')); 

load('S6_NO_FLOAT.mat');

%% EMG filtering
Trials=fieldnames(S6_NO_FLOAT);
Trials=Trials(9:11,1);

for Trialidx=1:length(Trials)
   TrialName=Trials{Trialidx}; 
   EMGSensors=fieldnames(S6_NO_FLOAT.(TrialName).Raw.EMG);
   if isstruct(S6_NO_FLOAT)
    for Sensoridx=1:length(EMGSensors)
        SensorName=EMGSensors{Sensoridx};
        if isstruct(S6_NO_FLOAT.(TrialName).Raw.EMG)
            %4-th order band pass filter (10-499 Hz)
            [bp1,bp2]=butter(4,[10 499]/(S6_NO_FLOAT.(TrialName).fsEMG/2)); 
            S6_NO_FLOAT.(TrialName).Filtered.(SensorName)=filtfilt(bp1,bp2,S6_NO_FLOAT.(TrialName).Raw.EMG.(SensorName));
        
            %rectification
            S6_NO_FLOAT.(TrialName).Rectified.(SensorName)=abs(S6_NO_FLOAT.(TrialName).Filtered.(SensorName));

            %4-th order Notch filter (50 Hz)
            [n1,n2]=butter(4,[49.8 50.2]/(S6_NO_FLOAT.(TrialName).fsEMG/2),'stop'); 
            S6_NO_FLOAT.(TrialName).Filtered3.(SensorName)=filtfilt(n1,n2,S6_NO_FLOAT.(TrialName).Rectified.(SensorName));
        
            %4-th order low-pass filter (10 Hz)
            [lp1,lp2]=butter(4,10/(S6_NO_FLOAT.(TrialName).fsEMG/2),'low'); 
            S6_NO_FLOAT.(TrialName).Filtered4.(SensorName)=filtfilt(lp1,lp2,S6_NO_FLOAT.(TrialName).Filtered3.(SensorName));
        end
    end
   end
end

%% Divide in gait cycles

%divide gait cycle 
%

gait_cycle1=200;
gait_cycle2=gait_cycle1+119;
gait_cycle3=gait_cycle2+129;
gait_cycle4=gait_cycle3+109;
vector1=S6_NO_FLOAT.T_01.Raw.Kin.LHIP-S6_NO_FLOAT.T_01.Raw.Kin.LKNE;


plot(vector1);
hold on 
scatter(100,vector1(100));
hold on 
scatter(gait_cycle1, vector1(gait_cycle1));
hold on 
scatter(gait_cycle2, vector1(gait_cycle2));
hold on
scatter(gait_cycle3, vector1(gait_cycle3));
hold on
scatter(gait_cycle4, vector1(gait_cycle4));

%vector2=
%vector3=




















