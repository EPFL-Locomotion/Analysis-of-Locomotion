clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Separating into gait cycles'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\GAIT FILES'));
addpath(fullfile('Separating into gait cycles')); 

load('NO_FLOAT_CRUTCHES.mat');

data1 = readtable('SCI_HCU_20150505_02OVGa_AD_01_GAIT.csv');
data2 = readtable('SCI_HCU_20150505_02OVGa_AD_02_GAIT.csv');
data3 = readtable('SCI_HCU_20150505_02OVGa_AD_03_GAIT.csv');

%% Getting the StepPoints

StepPoints.T_01.left=table2array(data1(1:8,4))*100;
StepPoints.T_01.right=table2array(data1(18:25,4))*100;
StepPoints.T_02.left=table2array(data2(1:8,4))*100;
StepPoints.T_02.right=table2array(data2(16:23,4))*100;
StepPoints.T_03.left=table2array(data3(1:8,4))*100;
StepPoints.T_03.right=table2array(data3(16:23,4))*100;

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


%% Eliminating some zeros values that will mess up everything

all_sensors={'LANK', 'LTOE', 'LWRA','LKNE','RANK', 'RTOE', 'RWRA', 'RKNE'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for i=1:length(all_sensors)
    
        NO_FLOAT_CRUTCHES.(trials{k}).Raw.Kin.(all_sensors{i})( all(~NO_FLOAT_CRUTCHES.(trials{k}).Raw.Kin.(all_sensors{i}),2), : ) = [];
       
    end
end

%% making struct for each cycles and sensor

sensors_left={'LANK', 'LTOE', 'LWRA','LKNE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(StepPoints.(trials{k}).left)-1
        for i=1:length(sensors_left)
            
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_left{i})=NO_FLOAT_CRUTCHES.(trials{k}).Raw.Kin.(sensors_left{i})(StepPoints.(trials{k}).left(n):StepPoints.(trials{k}).left(n+1),:);
        end
    end
end

sensors_right={'RANK', 'RTOE', 'RWRA', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(StepPoints.(trials{k}).right)-1
        for i=1:length(sensors_right)
            
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_right{i})=NO_FLOAT_CRUTCHES.(trials{k}).Raw.Kin.(sensors_right{i})(StepPoints.(trials{k}).right(n):StepPoints.(trials{k}).right(n+1),:);
            
        end
    end
end


%% Distance parameters

sensors={'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'};
trials={'T_01', 'T_02', 'T_03'};

% divide by fsKin to get it in seconds

for k=1:length(trials)
    for i=1:length(StepPoints.(trials{k}).right)-1
        for j=1:length(sensors)
            
            % stance
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.Stance.(sensors{j})=find(NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)==min(NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)))/100; 
            
            % swing
            A=NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j});
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.Swing.(sensors{j})=(find(A(:,3)==max(A(150:end,3)))-find(A(:,3)==min(A(100:end,3))))/100;

            % max toe step height
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.StepHeight.(sensors{j})=max(NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
       end
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(StepPoints.(trials{k}).right)-1
        
        for j=1:length(sensors)
         
            % max knee height
            NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.KneeHeight.(sensors{j})=max(NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
            
        end
    end
end

%% angles for both LEFT and RIGHT feet

sensors_left={'LTOE', 'LANK', 'LKNE', 'LWRA', 'LeftAngles'};
sensors_right={'RTOE', 'RANK', 'RKNE', 'RWRA', 'RightAngles'};
sensors=[sensors_left; sensors_right];

numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'};
trials={'T_01', 'T_02', 'T_03'};

for n=1:2
    
    for k=1:length(trials)
        
        for j=1:length(StepPoints.(trials{k}).left)-1
            
            for i=1:size(NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1}),1)
                
                TOE_points(i,:)=[NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,2) NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,3)];
                ANK_points(i,:)=[NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,2) NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,3)];
                KNE_points(i,:)=[NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,2) NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,3)];
                HIP_points(i,:)=[NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,2) NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,3)];
                
                VectorANKLE_TOE(i,:)=[(ANK_points(i,1)-TOE_points(i,1)) (ANK_points(i,2)-TOE_points(i,2))];
                VectorKNEE_ANKLE(i,:)=[(KNE_points(i,1)-ANK_points(i,1)) (KNE_points(i,2)-ANK_points(i,2))];
                VectorHIP_KNEE(i,:)=[(HIP_points(i,1)-KNE_points(i,1)) (HIP_points(i,2)-KNE_points(i,2))];
                
                Vertical = [0 1];
                Horizontal = [1 0];
                
                % hip angle to get extension/flexion
                Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).HipAngle(i)=Angle_hip_vertical(i)*180/pi;
                
                % joint angle knee
                Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeJointAngle(i)=Angle_knee(i)*180/pi;
                
                % elevation angle of the knee joint
                Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
                
                % joint angle ankle
                Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
                
                % elevation angle af the ankle joint
                Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
                
                % elevation toe from the ground
                Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
                NO_FLOAT_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
                
            end
        end
    end
end



