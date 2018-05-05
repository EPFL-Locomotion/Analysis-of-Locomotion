clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Separating into gait cycles'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\GAIT FILES'));
addpath(fullfile('Separating into gait cycles')); 

load('FLOAT_NO_CRUTCHES.mat');

%% EMG filtering
Trials=fieldnames(FLOAT_NO_CRUTCHES);
Trials=Trials(9:11,1);

for Trialidx=1:length(Trials)
   TrialName=Trials{Trialidx}; 
   EMGSensors=fieldnames(FLOAT_NO_CRUTCHES.(TrialName).Raw.EMG);
   if isstruct(FLOAT_NO_CRUTCHES)
    for Sensoridx=1:length(EMGSensors)
        SensorName=EMGSensors{Sensoridx};
        if isstruct(FLOAT_NO_CRUTCHES.(TrialName).Raw.EMG)
            %4-th order band pass filter (10-499 Hz)
            [bp1,bp2]=butter(4,[10 499]/(FLOAT_NO_CRUTCHES.(TrialName).fsEMG/2)); 
            FLOAT_NO_CRUTCHES.(TrialName).Filtered.(SensorName)=filtfilt(bp1,bp2,FLOAT_NO_CRUTCHES.(TrialName).Raw.EMG.(SensorName));
        
            %rectification
            FLOAT_NO_CRUTCHES.(TrialName).Rectified.(SensorName)=abs(FLOAT_NO_CRUTCHES.(TrialName).Filtered.(SensorName));

            %4-th order Notch filter (50 Hz)
            [n1,n2]=butter(4,[49.8 50.2]/(FLOAT_NO_CRUTCHES.(TrialName).fsEMG/2),'stop'); 
            FLOAT_NO_CRUTCHES.(TrialName).Filtered3.(SensorName)=filtfilt(n1,n2,FLOAT_NO_CRUTCHES.(TrialName).Rectified.(SensorName));
        
            %4-th order low-pass filter (10 Hz)
            [lp1,lp2]=butter(4,10/(FLOAT_NO_CRUTCHES.(TrialName).fsEMG/2),'low'); 
            FLOAT_NO_CRUTCHES.(TrialName).Filtered4.(SensorName)=filtfilt(lp1,lp2,FLOAT_NO_CRUTCHES.(TrialName).Filtered3.(SensorName));
        end
    end
   end
end



%% Divide in gait cycles





FLOAT_NO_CRUTCHES.T01.Gaitcycles.Cycle1.












