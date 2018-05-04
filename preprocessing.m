clear all
close all
clc

%Healthy without float

addpath(genpath('Healthy Recordings\Subject6\NO_FLOAT'));

load('S6_NO_FLOAT.mat');

Trials=fieldnames(S6_NO_FLOAT);
Trials=Trials(9:11,1);
EMGSensors=fieldnames(S6_NO_FLOAT.T_01.Raw.EMG);

for Trialidx=1:length(Trials)
   TrialName=Trials{Trialidx}; 
   if isstruct(S6_NO_FLOAT)
    for Sensoridx=1:length(EMGSensors)
        SensorName=EMGSensors{Sensoridx};
        if isstruct(S6_NO_FLOAT.(TrialName).Raw.EMG)
            %4-th order band pass filter (10-499 Hz)
            [bp1,bp2]=butter(2,[10 499]/(S6_NO_FLOAT.(TrialName).fsEMG/2)); 
            S6_NO_FLOAT.(TrialName).Filtered.(SensorName)=filter(bp1,bp2,S6_NO_FLOAT.(TrialName).Raw.EMG.(SensorName));
        
            %rectification
            S6_NO_FLOAT.(TrialName).Rectified.(SensorName)=abs(S6_NO_FLOAT.(TrialName).Filtered.(SensorName));

            %4-th order high-pass filter (30 Hz)
            [hp1,hp2]=butter(4,30/(S6_NO_FLOAT.(TrialName).fsEMG/2),'high'); 
            S6_NO_FLOAT.(TrialName).Filtered2.(SensorName)=filter(hp1,hp2,S6_NO_FLOAT.(TrialName).Rectified.(SensorName));
        
            %4-th order Notch filter (50 Hz)
            [n1,n2]=butter(4,[49.8 50.2]/(S6_NO_FLOAT.(TrialName).fsEMG/2),'stop'); 
            S6_NO_FLOAT.(TrialName).Filtered3.(SensorName)=filter(n1,n2,S6_NO_FLOAT.(TrialName).Filtered2.(SensorName));
        
            %4-th order low-pass filter (10 Hz)
            [lp1,lp2]=butter(4,10/(S6_NO_FLOAT.(TrialName).fsEMG/2),'low'); 
            S6_NO_FLOAT.(TrialName).Filtered4.(SensorName)=filter(lp1,lp2,S6_NO_FLOAT.(TrialName).Filtered3.(SensorName));
        end
    end
   end
end
subplot(3,2,1)
plot(S6_NO_FLOAT.T_01.Raw.EMG.LMG)
title('Raw signal');
subplot(3,2,2)
plot(S6_NO_FLOAT.T_01.Filtered.LMG)
title('BP filter 10-499 Hz');
subplot(3,2,3)
plot(S6_NO_FLOAT.T_01.Rectified.LMG)
title('Rectification');
subplot(3,2,4)
plot(S6_NO_FLOAT.T_01.Filtered2.LMG)
title('HP filter 30 Hz');
subplot(3,2,5)
plot(S6_NO_FLOAT.T_01.Filtered3.LMG)
title('Notch filter 50 Hz');
subplot(3,2,6)
plot(S6_NO_FLOAT.T_01.Filtered4.LMG)
title('LP filter 10 Hz');
