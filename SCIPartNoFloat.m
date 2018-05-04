clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Separating into gait cycles'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\GAIT FILES'));
addpath(fullfile('Separating into gait cycles')); 

load('NO_FLOAT_CRUTCHES.mat');

%% EMG filtering
Trials=fieldnames(NO_FLOAT_CRUTCHES);
Trials=Trials(9:11,1);

for Trialidx=1:length(Trials)
   TrialName=Trials{Trialidx}; 
   EMGSensors=fieldnames(NO_FLOAT_CRUTCHES.(TrialName).Raw.EMG);
   if isstruct(NO_FLOAT_CRUTCHES)
    for Sensoridx=1:length(EMGSensors)
        SensorName=EMGSensors{Sensoridx};
        if isstruct(NO_FLOAT_CRUTCHES.(TrialName).Raw.EMG)
            %4-th order band pass filter (10-499 Hz)
            [bp1,bp2]=butter(4,[10 499]/(NO_FLOAT_CRUTCHES.(TrialName).fsEMG/2)); 
            NO_FLOAT_CRUTCHES.(TrialName).Filtered.(SensorName)=filtfilt(bp1,bp2,NO_FLOAT_CRUTCHES.(TrialName).Raw.EMG.(SensorName));
        
            %rectification
            NO_FLOAT_CRUTCHES.(TrialName).Rectified.(SensorName)=abs(NO_FLOAT_CRUTCHES.(TrialName).Filtered.(SensorName));

            %4-th order Notch filter (50 Hz)
            [n1,n2]=butter(4,[49.8 50.2]/(NO_FLOAT_CRUTCHES.(TrialName).fsEMG/2),'stop'); 
            NO_FLOAT_CRUTCHES.(TrialName).Filtered3.(SensorName)=filtfilt(n1,n2,NO_FLOAT_CRUTCHES.(TrialName).Rectified.(SensorName));
        
            %4-th order low-pass filter (10 Hz)
            [lp1,lp2]=butter(4,10/(NO_FLOAT_CRUTCHES.(TrialName).fsEMG/2),'low'); 
            NO_FLOAT_CRUTCHES.(TrialName).Filtered4.(SensorName)=filtfilt(lp1,lp2,NO_FLOAT_CRUTCHES.(TrialName).Filtered3.(SensorName));
        end
    end
   end
end


%% Divide in gait cycles

