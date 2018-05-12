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

%% plot 3D left

curveLTOE=animatedline('Color', 'b');
curveLANK=animatedline('Color', 'r');
curveLHIP=animatedline('Color', 'm');
curveLKNE=animatedline('Color', 'g');

figure(1);
% set(gca, 'Xlim', [-50 280], 'YLim', [-3300 2800], 'ZLim', [0 1100]);
view(-95,14);
grid on;
hold on;
legend('LTOE', 'LANK', 'LHIP', 'LKNE');
xlabel('x');
ylabel('y');
zlabel('z');
title('3D plot of left sensors in space');

for i=1:length(S6_NO_FLOAT.T_01.Raw.Kin.LTOE)
    
    addpoints(curveLTOE, S6_NO_FLOAT.T_01.Raw.Kin.LTOE(i,1), S6_NO_FLOAT.T_01.Raw.Kin.LTOE(i,2), S6_NO_FLOAT.T_01.Raw.Kin.LTOE(i,3));
    addpoints(curveLANK, S6_NO_FLOAT.T_01.Raw.Kin.LANK(i,1), S6_NO_FLOAT.T_01.Raw.Kin.LANK(i,2), S6_NO_FLOAT.T_01.Raw.Kin.LANK(i,3));
    % plot3([S6_FLOAT.T_01.Raw.Kin.LTOE(i,1), S6_FLOAT.T_01.Raw.Kin.LANK(i,1)], [S6_FLOAT.T_01.Raw.Kin.LTOE(i,2), S6_FLOAT.T_01.Raw.Kin.LANK(i,2)], [S6_FLOAT.T_01.Raw.Kin.LTOE(i,3), S6_FLOAT.T_01.Raw.Kin.LANK(i,3)]);
   
    addpoints(curveLHIP, S6_NO_FLOAT.T_01.Raw.Kin.LHIP(i,1), S6_NO_FLOAT.T_01.Raw.Kin.LHIP(i,2), S6_NO_FLOAT.T_01.Raw.Kin.LHIP(i,3));
    addpoints(curveLKNE, S6_NO_FLOAT.T_01.Raw.Kin.LKNE(i,1), S6_NO_FLOAT.T_01.Raw.Kin.LKNE(i,2), S6_NO_FLOAT.T_01.Raw.Kin.LKNE(i,3));
   
    drawnow limitrate
    % pause(0.1);
    
end
hold off;

%% Eliminating some zeros values that will mess up everything

all_sensors={'LANK', 'LTOE', 'LHIP','LKNE','RANK', 'RTOE', 'RHIP', 'RKNE'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for i=1:length(all_sensors)
    
        S6_NO_FLOAT.(trials{k}).Raw.Kin.(all_sensors{i})( all(~S6_NO_FLOAT.(trials{k}).Raw.Kin.(all_sensors{i}),2), : ) = [];
        
        
    end
end

%% split in cycles right T_01

IC_right(1)=max(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(1:150,3));
IC_right(2)=max(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(151:220,3));
IC_right(3)=max(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(221:400,3));
IC_right(4)=max(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(401:500,3));
IC_right(5)=max(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(501:end,3));

for i=1:5
    
PosIC_right.T_01(i)=find(S6_NO_FLOAT.T_01.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% split in cycles left T_01

IC_left(1)=max(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(1:200,3));
IC_left(2)=max(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(201:320,3));
IC_left(3)=max(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(321:450,3));
IC_left(4)=max(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(451:550,3));
IC_left(5)=max(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(551:end,3));

for i=1:5
    
PosIC_left.T_01(i)=find(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles left T_02

IC_left(1)=max(S6_NO_FLOAT.T_02.Raw.Kin.LTOE(1:150,3));
IC_left(2)=max(S6_NO_FLOAT.T_02.Raw.Kin.LTOE(151:300,3));
IC_left(3)=max(S6_NO_FLOAT.T_02.Raw.Kin.LTOE(301:400,3));
IC_left(4)=max(S6_NO_FLOAT.T_02.Raw.Kin.LTOE(401:500,3));

for i=1:4
    
PosIC_left.T_02(i)=find(S6_NO_FLOAT.T_02.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles right T_02

IC_right(1)=max(S6_NO_FLOAT.T_02.Raw.Kin.RTOE(1:200,3));
IC_right(2)=max(S6_NO_FLOAT.T_02.Raw.Kin.RTOE(201:300,3));
IC_right(3)=max(S6_NO_FLOAT.T_02.Raw.Kin.RTOE(301:400,3));
IC_right(4)=max(S6_NO_FLOAT.T_02.Raw.Kin.RTOE(401:end,3));

for i=1:4
    
PosIC_right.T_02(i)=find(S6_NO_FLOAT.T_02.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% split in cycles left T_03

IC_left(1)=max(S6_NO_FLOAT.T_03.Raw.Kin.LTOE(1:150,3));
IC_left(2)=max(S6_NO_FLOAT.T_03.Raw.Kin.LTOE(151:250,3));
IC_left(3)=max(S6_NO_FLOAT.T_03.Raw.Kin.LTOE(251:400,3));
IC_left(4)=max(S6_NO_FLOAT.T_03.Raw.Kin.LTOE(401:500,3));

for i=1:4
    
PosIC_left.T_03(i)=find(S6_NO_FLOAT.T_03.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles right T_03

IC_right(1)=max(S6_NO_FLOAT.T_03.Raw.Kin.RTOE(1:200,3));
IC_right(2)=max(S6_NO_FLOAT.T_03.Raw.Kin.RTOE(201:300,3));
IC_right(3)=max(S6_NO_FLOAT.T_03.Raw.Kin.RTOE(301:400,3));
IC_right(4)=max(S6_NO_FLOAT.T_03.Raw.Kin.RTOE(401:500,3));

for i=1:4
    
PosIC_right.T_03(i)=find(S6_NO_FLOAT.T_03.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% making struct for each cycles and sensor

sensors_left={'LANK', 'LTOE', 'LHIP','LKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(PosIC_left.(trials{k}))-1
        for i=1:length(sensors_left)
            
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_left{i})=S6_NO_FLOAT.(trials{k}).Raw.Kin.(sensors_left{i})(PosIC_left.(trials{k})(n):PosIC_left.(trials{k})(n+1),:);
        end
    end
end

sensors_right={'RANK', 'RTOE', 'RHIP', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(PosIC_left.(trials{k}))-1
        for i=1:length(sensors_right)
            
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_right{i})=S6_NO_FLOAT.(trials{k}).Raw.Kin.(sensors_right{i})(PosIC_right.(trials{k})(n):PosIC_right.(trials{k})(n+1),:);
            
        end
    end
end

%% Distance parameters

sensors={'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

% divide by fsKin to get it in seconds

for k=1:length(trials)
    for i=1:length(PosIC_left.(trials{k}))-1
        for j=1:length(sensors)
            
            % stance
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.Stance.(sensors{j})=find(S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)==min(S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)))/100; 
            
            % swing
            A=S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j});
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.Swing.(sensors{j})=(find(A(:,3)==max(A(60:end,3)))-find(A(:,3)==min(A(50:end,3))))/100;

            % max toe step height
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.StepHeight.(sensors{j})=max(S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
       end
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(PosIC_left.(trials{k}))-1
        
        for j=1:length(sensors)
         
            % max knee height
            S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.KneeHeight.(sensors{j})=max(S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
            
        end
    end
end

%% angles for both LEFT and RIGHT feet

sensors_left={'LTOE', 'LANK', 'LKNE', 'LHIP', 'LeftAngles'};
sensors_right={'RTOE', 'RANK', 'RKNE', 'RHIP', 'RightAngles'};
sensors=[sensors_left; sensors_right];

numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for n=1:2
    
    for k=1:length(trials)
        
        for j=1:length(PosIC_left.(trials{k}))-1
            
            for i=1:size(S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1}),1)
                
                TOE_points(i,:)=[S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,2) S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,3)];
                ANK_points(i,:)=[S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,2) S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,3)];
                KNE_points(i,:)=[S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,2) S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,3)];
                HIP_points(i,:)=[S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,2) S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,3)];
                
                VectorANKLE_TOE(i,:)=[(ANK_points(i,1)-TOE_points(i,1)) (ANK_points(i,2)-TOE_points(i,2))];
                VectorKNEE_ANKLE(i,:)=[(KNE_points(i,1)-ANK_points(i,1)) (KNE_points(i,2)-ANK_points(i,2))];
                VectorHIP_KNEE(i,:)=[(HIP_points(i,1)-KNE_points(i,1)) (HIP_points(i,2)-KNE_points(i,2))];
                
                Vertical = [0 1];
                Horizontal = [1 0];
                
                % hip angle to get extension/flexion
                Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).HipAngle(i)=Angle_hip_vertical(i)*180/pi;
                
                % joint angle knee
                Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeJointAngle(i)=Angle_knee(i)*180/pi;
                
                % elevation angle of the knee joint
                Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
                
                % joint angle ankle
                Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
                
                % elevation angle af the ankle joint
                Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
                
                % elevation toe from the ground
                Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
                
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
     for j=1:length(PosIC_left.(trials{k}))-1
        for l=1:length(EMGType)   
            for i=1:length(EMGSensorsLeft)
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsLeft{i})=S6_NO_FLOAT.(trials{k}).(EMGType{l}).(EMGSensorsLeft{i})((10*PosIC_left.(trials{k})(j)):(10*PosIC_left.(trials{k})(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(PosIC_left.(trials{k}))-1
         for i=1:length(EMGSensorsLeft)
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsLeft{i})=S6_NO_FLOAT.(trials{k}).Raw.EMG.(EMGSensorsLeft{i})((10*PosIC_left.(trials{k})(j)):(10*PosIC_left.(trials{k})(j+1)-1));
         end
     end
end




%divide into gait cycles - right
for k=1:length(trials)
     for j=1:length(PosIC_right.(trials{k}))-1
        for l=1:length(EMGType)   
            for i=1:length(EMGSensorsRight)
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsRight{i})=S6_NO_FLOAT.(trials{k}).(EMGType{l}).(EMGSensorsRight{i})((10*PosIC_right.(trials{k})(j)):(10*PosIC_right.(trials{k})(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(PosIC_right.(trials{k}))-1
         for i=1:length(EMGSensorsRight)
                S6_NO_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsRight{i})=S6_NO_FLOAT.(trials{k}).Raw.EMG.(EMGSensorsRight{i})((10*PosIC_right.(trials{k})(j)):(10*PosIC_right.(trials{k})(j+1)-1));
         end
     end
end


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_NO_FLOAT.(trials{k}).Raw.EMG.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Raw signal');


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_NO_FLOAT.(trials{k}).Filtered.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Filtered signal');


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_NO_FLOAT.(trials{k}).Filtered4.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Final filtered signal');


%--
figure
     for j=2:3
         for i=1:length(fieldnames(S6_NO_FLOAT.(trials{3}).GaitCycles))
            subplot(4,4,4*(j-1)+i)
            plot(S6_NO_FLOAT.(trials{3}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
         end   
     end
suptitle('Trial 3 - Raw');


figure
     for j=2:3
         for i=1:length(fieldnames(S6_NO_FLOAT.(trials{1}).GaitCycles))
            subplot(4,4,4*(j-1)+i)
            plot(S6_NO_FLOAT.(trials{1}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
         end   
     end
suptitle('Trial 1 - Final filtered');


%%Bursts calculation visually
Bursts.Position.T_01.GaitCycles.One.LTA=[108 705];
Bursts.Position.T_01.GaitCycles.Two.LTA=[80 750];
Bursts.Position.T_01.GaitCycles.Three.LTA=[44 736];
Bursts.Position.T_01.GaitCycles.Four.LTA=[90 761];

Bursts.Position.T_01.GaitCycles.One.RMG=[80 897];
Bursts.Position.T_01.GaitCycles.Two.RMG=[72 752];
Bursts.Position.T_01.GaitCycles.Three.RMG=[150 750];
Bursts.Position.T_01.GaitCycles.Four.RMG=[147 824];

Bursts.Position.T_02.GaitCycles.One.LTA=[193 897];
Bursts.Position.T_02.GaitCycles.Two.LTA=[146 700];
Bursts.Position.T_02.GaitCycles.Three.LTA=[110 710];

Bursts.Position.T_02.GaitCycles.One.RMG=[86 750];
Bursts.Position.T_02.GaitCycles.Two.RMG=[100 711];
Bursts.Position.T_02.GaitCycles.Three.RMG=[94 700];


Bursts.Position.T_03.GaitCycles.One.LTA=[390 1805];
Bursts.Position.T_03.GaitCycles.Two.LTA=[260 1773];
Bursts.Position.T_03.GaitCycles.Three.LTA=[260 1773];

Bursts.Position.T_03.GaitCycles.One.RMG=[535 1900];
Bursts.Position.T_03.GaitCycles.Two.RMG=[458 1805];
Bursts.Position.T_03.GaitCycles.Three.RMG=[260 1773];

%%Calculation of EMG parameters
for i=1:3
    for j=1:length(fieldnames(S6_NO_FLOAT.(trials{i}).GaitCycles))
        for k=2:3
            %Each burst falls between two gait cycles, i.e. each gait
            %cycles has two half-bursts (one at the beginning and one at
            %the end)
            end1burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(1);%end of the first half-burst
            start2burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(2);%start of the second half-burst
            lengthburst=length(S6_NO_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensors{k}));
        
            %Duration in seconds (i.e. divided by 1000), i.e. how much of the gait cycle is burst
            Bursts.Duration.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=0.001*(end1burst-1+lengthburst-start2burst);
        
            %Max value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MaxValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=max(S6_NO_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k}));
            
            %Mean value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MeanValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=mean(S6_NO_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k})([1:end1burst,start2burst:lengthburst]));
        end
    end
end















