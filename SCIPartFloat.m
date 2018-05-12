clear all
close all
clc

% adding the paths and loading data

addpath(genpath('SCI subject (updated file)'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\GAIT FILES'));
 
load('FLOAT_NO_CRUTCHES.mat');

data1 = readtable('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_03_GAIT.csv');
data2 = readtable('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_04_GAIT.csv');
data3 = readtable('SCI_HCU_20150505_04OVGb_45BWS_vFWD_noAD_05_GAIT.csv');

%% Getting the StepPoints

StepPoints.T_01.left=table2array(data1(1:6,4))*100;
StepPoints.T_01.right=table2array(data1(14:19,4))*100;
StepPoints.T_02.left=table2array(data2(1:5,4))*100;
StepPoints.T_02.right=table2array(data2(12:16,4))*100;
StepPoints.T_03.left=table2array(data3(1:6,4))*100;
StepPoints.T_03.right=table2array(data3(12:17,4))*100;

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



%% Eliminating some zeros values that will mess up everything

all_sensors={'LANK', 'LTOE', 'LWRA','LKNE','RANK', 'RTOE', 'RWRA', 'RKNE'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for i=1:length(all_sensors)
    
        FLOAT_NO_CRUTCHES.(trials{k}).Raw.Kin.(all_sensors{i})( all(~FLOAT_NO_CRUTCHES.(trials{k}).Raw.Kin.(all_sensors{i}),2), : ) = [];
       
    end
end


%% making struct for each cycles and sensor

sensors_left={'LANK', 'LTOE', 'LWRA','LKNE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(StepPoints.(trials{k}).left)-1
        for i=1:length(sensors_left)
            
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_left{i})=FLOAT_NO_CRUTCHES.(trials{k}).Raw.Kin.(sensors_left{i})(StepPoints.(trials{k}).left(n):StepPoints.(trials{k}).left(n+1),:);
        end
    end
end

sensors_right={'RANK', 'RTOE', 'RWRA', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(StepPoints.(trials{k}).right)-1
        for i=1:length(sensors_right)
            
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_right{i})=FLOAT_NO_CRUTCHES.(trials{k}).Raw.Kin.(sensors_right{i})(StepPoints.(trials{k}).right(n):StepPoints.(trials{k}).right(n+1),:);
            
        end
    end
end

%% Distance parameters

sensors={'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six'};
trials={'T_01', 'T_02', 'T_03'};

% divide by fsKin to get it in seconds

for k=1:length(trials)
    for i=1:length(StepPoints.(trials{k}).right)-1
        for j=1:length(sensors)
            
            % stance
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.Stance.(sensors{j})=find(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)==min(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)))/100; 
            
            % swing
            A=FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j});
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.Swing.(sensors{j})=(find(A(:,3)==max(A(150:end,3)))-find(A(:,3)==min(A(100:end,3))))/100;

            % max toe step height
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.StepHeight.(sensors{j})=max(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
       end
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(StepPoints.(trials{k}).right)-1
        
        for j=1:length(sensors)
         
            % max knee height
            FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.KneeHeight.(sensors{j})=max(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
            
        end
    end
end


%% angles for both LEFT and RIGHT feet

sensors_left={'LTOE', 'LANK', 'LKNE', 'LWRA', 'LeftAngles'};
sensors_right={'RTOE', 'RANK', 'RKNE', 'RWRA', 'RightAngles'};
sensors=[sensors_left; sensors_right];

numbers={'One', 'Two', 'Three', 'Four', 'Five', 'Six'};
trials={'T_01', 'T_02', 'T_03'};

for n=1:2
    
    for k=1:length(trials)
        
        for j=1:length(StepPoints.(trials{k}).left)-1
            
            for i=1:size(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1}),1)
                
                TOE_points(i,:)=[FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,2) FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,3)];
                ANK_points(i,:)=[FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,2) FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,3)];
                KNE_points(i,:)=[FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,2) FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,3)];
                HIP_points(i,:)=[FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,2) FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,3)];
                
                VectorANKLE_TOE(i,:)=[(ANK_points(i,1)-TOE_points(i,1)) (ANK_points(i,2)-TOE_points(i,2))];
                VectorKNEE_ANKLE(i,:)=[(KNE_points(i,1)-ANK_points(i,1)) (KNE_points(i,2)-ANK_points(i,2))];
                VectorHIP_KNEE(i,:)=[(HIP_points(i,1)-KNE_points(i,1)) (HIP_points(i,2)-KNE_points(i,2))];
                
                Vertical = [0 1];
                Horizontal = [1 0];
                
                % hip angle to get extension/flexion
                Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).HipAngle(i)=Angle_hip_vertical(i)*180/pi;
                
                % joint angle knee
                Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeJointAngle(i)=Angle_knee(i)*180/pi;
                
                % elevation angle of the knee joint
                Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
                
                % joint angle ankle
                Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
                
                % elevation angle af the ankle joint
                Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
                
                % elevation toe from the ground
                Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
                
            end
        end
    end
end



%% EMG parameters
EMGType={'Filtered','Rectified','Filtered3','Filtered4'};
EMGSensorsLeft={'LMG','LTA'};
EMGSensorsRight={'RMG','RTA'};


%divide into gait cycles - left
for k=1:length(trials)
     for j=1:length(StepPoints.(trials{k}).left)-1
        for l=1:length(EMGType)   
            for i=1:length(EMGSensorsLeft)
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsLeft{i})=FLOAT_NO_CRUTCHES.(trials{k}).(EMGType{l}).(EMGSensorsLeft{i})((10*StepPoints.(trials{k}).left(j)):(10*StepPoints.(trials{k}).left(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(StepPoints.(trials{k}).left)-1
         for i=1:length(EMGSensorsLeft)
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsLeft{i})=FLOAT_NO_CRUTCHES.(trials{k}).Raw.EMG.(EMGSensorsLeft{i})((10*StepPoints.(trials{k}).left(j)):(10*StepPoints.(trials{k}).left(j+1)-1));
         end
     end
end




%divide into gait cycles - right
for k=1:length(trials)
     for j=1:length(StepPoints.(trials{k}).right)-1
        for l=1:length(EMGType)   
            for i=1:length(EMGSensorsRight)
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsRight{i})=FLOAT_NO_CRUTCHES.(trials{k}).(EMGType{l}).(EMGSensorsRight{i})((10*StepPoints.(trials{k}).right(j)):(10*StepPoints.(trials{k}).right(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(StepPoints.(trials{k}).right)-1
         for i=1:length(EMGSensorsRight)
                FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsRight{i})=FLOAT_NO_CRUTCHES.(trials{k}).Raw.EMG.(EMGSensorsRight{i})((10*StepPoints.(trials{k}).right(j)):(10*StepPoints.(trials{k}).right(j+1)-1));
         end
     end
end


figure
for k=1:length(trials)
    cont=1;
     for i=[31,32,39,40]
         subplot(3,4,4*(k-1)+cont)
         plot(FLOAT_NO_CRUTCHES.(trials{k}).Raw.EMG.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
         cont=cont+1;
     end
end
suptitle('Raw signal');


figure
for k=1:length(trials)
    cont=1; 
    for i=[31,32,39,40]
         subplot(3,4,4*(k-1)+cont)
         plot(FLOAT_NO_CRUTCHES.(trials{k}).Filtered.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
         cont=cont+1;
    end
end
suptitle('Filtered signal');


figure
for k=1:length(trials)
    cont=1;
    for i=[31,32,39,40]
        subplot(3,4,4*(k-1)+cont)
        plot(FLOAT_NO_CRUTCHES.(trials{k}).Filtered4.(EMGSensors{i}));
        title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
        cont=cont+1;
    end
end
suptitle('Final filtered signal');


%--
figure
     for i=1:length(fieldnames(FLOAT_NO_CRUTCHES.(trials{3}).GaitCycles))
     cont=1;
         for j=[32,39]
            subplot(2,5,5*(cont-1)+i)
            plot(FLOAT_NO_CRUTCHES.(trials{3}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
            cont=cont+1;
         end
      end
suptitle('Trial 3 - Raw');


figure
     for i=1:length(fieldnames(FLOAT_NO_CRUTCHES.(trials{3}).GaitCycles))
     cont=1;
         for j=[32,39]
            subplot(2,5,5*(cont-1)+i)
            plot(FLOAT_NO_CRUTCHES.(trials{3}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
            cont=cont+1;
         end   
     end
suptitle('Trial 3 - Final filtered');


%% BURSTS CALCULATION FROM STD
%---LTA
for k=1:length(trials)
    for i=1:length(fieldnames(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles))
         for j=32%[32,39]
            %Assumption#1.1: for LTA between 500 and 600 it is always noise
            %Assumption#1.2: we put a higher and a lower limit for 2*std
            clear StdNoise;
            clear indeces;
            StdNoise=min([std(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j})(500:600)),0.0375/2]);
            StdNoise=max([StdNoise,0.01/2]);
            indeces=find(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j})>=(2*StdNoise));
            %Assumption#2: minimum burst length is 300
            Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=[];
            for l=1:(length(indeces)-1)
               if (indeces(l+1)-indeces(l)>300)  
                  Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=cat(1,Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}),indeces(l),indeces(l+1)); 
               end
            end
            
            %Particular case#1.a: no initial (i.e. let's say before first 25% of gait cycle) burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)>0.25*length((FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=1;
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=indeces(1);
            end    
            %Particular case#1.b: no final burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)<0.30*length((FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                 Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=indeces(length(indeces));
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=length(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
            end    
            
            %Particular case#2: false internal burst
            for l=1:(length(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}))-2)
               if (Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1)-Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l)<=300)  
                    %remove the two successive burst indeces (because that
                    %burst is just inside the real big burst)
                    Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}) = Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})~=Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1));
                    Bursts.Position.(trials{1}).GaitCycles.(numbers{i}).(EMGSensors{j}) = Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})~=Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1));
               end
            end
            
         end   
    end
end

%---RMG
for k=1:length(trials)
    for i=1:length(fieldnames(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles))
         for j=39%[32,39]
            %Assumption#1.1: for LTA between 1750 and 1850 it is always noise
            %Assumption#1.2: we put a higher and a lower limit for 2*std
            clear StdNoise;
            clear indeces;
            StdNoise=min([std(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j})(1750:1850)),0.0375/2]);
            indeces=find(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j})>=(2*StdNoise));
            %Assumption#2: minimum burst length is 300
            Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=[];
            for l=1:(length(indeces)-1)
               if (indeces(l+1)-indeces(l)>300)  
                  Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=cat(1,Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}),indeces(l),indeces(l+1)); 
               end
            end
            
            
            %Particular case#1.a: no initial (i.e. let's say before first 25% of gait cycle) burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)>0.25*length((FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=1;
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=indeces(1);
            end    
            %Particular case#1.b: no final burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)<0.30*length((FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=indeces(length(indeces));
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=length(FLOAT_NO_CRUTCHES.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
            end    
                        
            %Particular case#2: false internal burst
            for l=1:(length(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}))-2)
               if (Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1)-Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l)<=300)  
                    %remove the two successive burst indeces (because that
                    %burst is just inside the real big burst)
                    Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}) = Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})~=Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1));
                    Bursts.Position.(trials{1}).GaitCycles.(numbers{i}).(EMGSensors{j}) = Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})~=Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(l+1));
               end
            end
            
         end   
    end
end




%% 

% %%Bursts calculation visually
% 
% Bursts.Position.T_01.GaitCycles.One.LTA=[108 705];
% Bursts.Position.T_01.GaitCycles.Two.LTA=[80 750];
% Bursts.Position.T_01.GaitCycles.Three.LTA=[44 736];
% Bursts.Position.T_01.GaitCycles.Four.LTA=[90 761];
% 
% Bursts.Position.T_01.GaitCycles.One.RMG=[80 897];
% Bursts.Position.T_01.GaitCycles.Two.RMG=[72 752];
% Bursts.Position.T_01.GaitCycles.Three.RMG=[150 750];
% Bursts.Position.T_01.GaitCycles.Four.RMG=[147 824];
% 
% Bursts.Position.T_02.GaitCycles.One.LTA=[193 897];
% Bursts.Position.T_02.GaitCycles.Two.LTA=[146 700];
% Bursts.Position.T_02.GaitCycles.Three.LTA=[110 710];
% 
% Bursts.Position.T_02.GaitCycles.One.RMG=[86 750];
% Bursts.Position.T_02.GaitCycles.Two.RMG=[100 711];
% Bursts.Position.T_02.GaitCycles.Three.RMG=[94 700];
% 
% 
% Bursts.Position.T_03.GaitCycles.One.LTA=[390 1805];
% Bursts.Position.T_03.GaitCycles.Two.LTA=[260 1773];
% Bursts.Position.T_03.GaitCycles.Three.LTA=[260 1773];
% 
% Bursts.Position.T_03.GaitCycles.One.RMG=[535 1900];
% Bursts.Position.T_03.GaitCycles.Two.RMG=[458 1805];
% Bursts.Position.T_03.GaitCycles.Three.RMG=[260 1773];
% 
%%Calculation of EMG parameters
for i=1:3
    for j=1:length(fieldnames(FLOAT_NO_CRUTCHES.(trials{i}).GaitCycles))
        for k=[32,39]
            %Each burst falls between two gait cycles, i.e. each gait
            %cycles has two half-bursts (one at the beginning and one at
            %the end)
            end1burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(1);%end of the first half-burst
            start2burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(2);%start of the second half-burst
            lengthburst=length(FLOAT_NO_CRUTCHES.(trials{i}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensors{k}));
        
            %Duration in seconds (i.e. divided by 1000), i.e. how much of the gait cycle is burst
            Bursts.Duration.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=0.001*(end1burst-1+lengthburst-start2burst);
        
            %Max value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MaxValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=max(FLOAT_NO_CRUTCHES.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k}));
            
            %Mean value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MeanValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=mean(FLOAT_NO_CRUTCHES.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k})([1:end1burst,start2burst:lengthburst]));
        end
    end
end




