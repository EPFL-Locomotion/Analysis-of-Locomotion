clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Healthy Recordings\Subject6\FLOAT'));

load('S6_FLOAT.mat');

%% EMG filtering
Trials=fieldnames(S6_FLOAT);
Trials=Trials(9:11,1);

for Trialidx=1:length(Trials)
    TrialName=Trials{Trialidx};
    EMGSensors=fieldnames(S6_FLOAT.(TrialName).Raw.EMG);
    if isstruct(S6_FLOAT)
        for Sensoridx=1:length(EMGSensors)
            SensorName=EMGSensors{Sensoridx};
            if isstruct(S6_FLOAT.(TrialName).Raw.EMG)
                %4-th order band pass filter (10-499 Hz)
                [bp1,bp2]=butter(4,[10 499]/(S6_FLOAT.(TrialName).fsEMG/2));
                S6_FLOAT.(TrialName).Filtered.(SensorName)=filtfilt(bp1,bp2,S6_FLOAT.(TrialName).Raw.EMG.(SensorName));
                
                %rectification
                S6_FLOAT.(TrialName).Rectified.(SensorName)=abs(S6_FLOAT.(TrialName).Filtered.(SensorName));
                
                %4-th order Notch filter (50 Hz)
                [n1,n2]=butter(4,[49.8 50.2]/(S6_FLOAT.(TrialName).fsEMG/2),'stop');
                S6_FLOAT.(TrialName).Filtered3.(SensorName)=filtfilt(n1,n2,S6_FLOAT.(TrialName).Rectified.(SensorName));
                
                %4-th order low-pass filter (10 Hz)
                [lp1,lp2]=butter(4,10/(S6_FLOAT.(TrialName).fsEMG/2),'low');
                S6_FLOAT.(TrialName).Filtered4.(SensorName)=filtfilt(lp1,lp2,S6_FLOAT.(TrialName).Filtered3.(SensorName));
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
set(gca, 'Xlim', [-50 280], 'YLim', [-3300 2800], 'ZLim', [0 1100]);
view(-95,14);
grid on;
hold on;
legend('LTOE', 'LANK', 'LHIP', 'LKNE');
xlabel('x');
ylabel('y');
zlabel('Z');
title('3D plot of left sensors in space');

for i=1:1033
    
    addpoints(curveLTOE, S6_FLOAT.T_01.Raw.Kin.LTOE(i,1), S6_FLOAT.T_01.Raw.Kin.LTOE(i,2), S6_FLOAT.T_01.Raw.Kin.LTOE(i,3));
    addpoints(curveLANK, S6_FLOAT.T_01.Raw.Kin.LANK(i,1), S6_FLOAT.T_01.Raw.Kin.LANK(i,2), S6_FLOAT.T_01.Raw.Kin.LANK(i,3));
    % plot3([S6_FLOAT.T_01.Raw.Kin.LTOE(i,1), S6_FLOAT.T_01.Raw.Kin.LANK(i,1)], [S6_FLOAT.T_01.Raw.Kin.LTOE(i,2), S6_FLOAT.T_01.Raw.Kin.LANK(i,2)], [S6_FLOAT.T_01.Raw.Kin.LTOE(i,3), S6_FLOAT.T_01.Raw.Kin.LANK(i,3)]);
   
    addpoints(curveLHIP, S6_FLOAT.T_01.Raw.Kin.LHIP(i,1), S6_FLOAT.T_01.Raw.Kin.LHIP(i,2), S6_FLOAT.T_01.Raw.Kin.LHIP(i,3));
    addpoints(curveLKNE, S6_FLOAT.T_01.Raw.Kin.LKNE(i,1), S6_FLOAT.T_01.Raw.Kin.LKNE(i,2), S6_FLOAT.T_01.Raw.Kin.LKNE(i,3));
   
    drawnow limitrate
    % pause(0.1);
    
end
hold off;

%% Eliminating some zeros values that will mess up everything

all_sensors={'LANK', 'LTOE', 'LHIP','LKNE','RANK', 'RTOE', 'RHIP', 'RKNE'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for i=1:length(all_sensors)
    
        S6_FLOAT.(trials{k}).Raw.Kin.(all_sensors{i})( all(~S6_FLOAT.(trials{k}).Raw.Kin.(all_sensors{i}),2), : ) = [];
        
        
    end
end

%             S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).(sensors_left{i})=S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).(sensors_left{i})(find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).(sensors_left{i})~=0));
%            S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).(sensors_left{i})( all(~S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).(sensors_left{i}),2), : ) = [];
        

%% split in cycles right T_01

IC_right(1)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(1:250,3));
IC_right(2)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(251:500,3));
IC_right(3)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(501:700,3));
IC_right(4)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(701:900,3));
IC_right(5)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(901:end,3));

for i=1:5
    
PosIC_right.T_01(i)=find(S6_FLOAT.T_01.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% split in cycles left T_01

IC_left(1)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(1:50,3));
IC_left(2)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(51:400,3));
IC_left(3)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(401:600,3));
IC_left(4)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(601:800,3));
IC_left(5)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(800:end,3));


for i=1:5
    
PosIC_left.T_01(i)=find(S6_FLOAT.T_01.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles left T_02

IC_left(1)=max(S6_FLOAT.T_02.Raw.Kin.LTOE(1:200,3));
IC_left(2)=max(S6_FLOAT.T_02.Raw.Kin.LTOE(201:500,3));
IC_left(3)=max(S6_FLOAT.T_02.Raw.Kin.LTOE(501:700,3));
IC_left(4)=max(S6_FLOAT.T_02.Raw.Kin.LTOE(701:end,3));

for i=1:4
    
PosIC_left.T_02(i)=find(S6_FLOAT.T_02.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles right T_02

IC_right(1)=max(S6_FLOAT.T_02.Raw.Kin.RTOE(1:200,3));
IC_right(2)=max(S6_FLOAT.T_02.Raw.Kin.RTOE(201:470,3));
IC_right(3)=max(S6_FLOAT.T_02.Raw.Kin.RTOE(471:700,3));
IC_right(4)=max(S6_FLOAT.T_02.Raw.Kin.RTOE(701:end,3));

for i=1:4
    
PosIC_right.T_02(i)=find(S6_FLOAT.T_02.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% split in cycles left T_03

IC_left(1)=max(S6_FLOAT.T_03.Raw.Kin.LTOE(50:300,3));
IC_left(2)=max(S6_FLOAT.T_03.Raw.Kin.LTOE(301:550,3));
IC_left(3)=max(S6_FLOAT.T_03.Raw.Kin.LTOE(551:end,3));

for i=1:3
    
PosIC_left.T_03(i)=find(S6_FLOAT.T_03.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% split in cycles right T_03

IC_right(1)=max(S6_FLOAT.T_03.Raw.Kin.RTOE(1:200,3));
IC_right(2)=max(S6_FLOAT.T_03.Raw.Kin.RTOE(201:450,3));
IC_right(3)=max(S6_FLOAT.T_03.Raw.Kin.RTOE(451:end,3));

for i=1:3
    
PosIC_right.T_03(i)=find(S6_FLOAT.T_03.Raw.Kin.RTOE(:,3)==IC_right(i));

end


%% making struct for each cycles and sensor

sensors_left={'LANK', 'LTOE', 'LHIP','LKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(PosIC_left.(trials{k}))-1
        for i=1:length(sensors_left)
            
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_left{i})=S6_FLOAT.(trials{k}).Raw.Kin.(sensors_left{i})(PosIC_left.(trials{k})(n):PosIC_left.(trials{k})(n+1),:);
        end
    end
end

sensors_right={'RANK', 'RTOE', 'RHIP', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for n=1:length(PosIC_left.(trials{k}))-1
        for i=1:length(sensors_right)
            
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{n}).Kin.(sensors_right{i})=S6_FLOAT.(trials{k}).Raw.Kin.(sensors_right{i})(PosIC_right.(trials{k})(n):PosIC_right.(trials{k})(n+1),:);
            
        end
    end
end



%% sample plots 
%(avoiding to plot the x direction cause we suppose he moves straigh)

figure;
plot(S6_FLOAT.T_01.GaitCycles.One.Kin.LANK(:,2), S6_FLOAT.T_01.GaitCycles.One.Kin.LANK(:,3)); 

figure;
plot(S6_FLOAT.T_01.GaitCycles.Two.Kin.RTOE(:,2), S6_FLOAT.T_01.GaitCycles.Two.Kin.RTOE(:,3)); 



%% Distance parameters

sensors={'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

% divide by fsKin to get it in seconds

for k=1:length(trials)
    for i=1:length(PosIC_left.(trials{k}))-1
        for j=1:length(sensors)
            
            % stance
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.Stance.(sensors{j})=find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)==min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3)))/100; 
            
            % swing
            A=S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j});
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.Swing.(sensors{j})=(find(A(:,3)==max(A(150:end,3)))-find(A(:,3)==min(A(100:end,3))))/100;

            % max toe step height
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.StepHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
       end
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(PosIC_left.(trials{k}))-1
        
        for j=1:length(sensors)
         
            % max knee height
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.KneeHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Kin.(sensors{j})(:,3));
            
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
            
            for i=1:size(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1}),1)
                
                TOE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,1})(i,3)];
                ANK_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,2})(i,3)];
                KNE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,3})(i,3)];
                HIP_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,4})(i,3)];
                
                VectorANKLE_TOE(i,:)=[(ANK_points(i,1)-TOE_points(i,1)) (ANK_points(i,2)-TOE_points(i,2))];
                VectorKNEE_ANKLE(i,:)=[(KNE_points(i,1)-ANK_points(i,1)) (KNE_points(i,2)-ANK_points(i,2))];
                VectorHIP_KNEE(i,:)=[(HIP_points(i,1)-KNE_points(i,1)) (HIP_points(i,2)-KNE_points(i,2))];
                
                Vertical = [0 1];
                Horizontal = [1 0];
                
                % hip angle to get extension/flexion
                Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).HipAngle(i)=Angle_hip_vertical(i)*180/pi;
                
                % joint angle knee
                Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeJointAngle(i)=Angle_knee(i)*180/pi;
                
                % elevation angle of the knee joint
                Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
                
                % joint angle ankle
                Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
                
                % elevation angle af the ankle joint
                Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
                
                % elevation toe from the ground
                Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.(sensors{n,5}).ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
                
            end
        end
    end
end



%% ALL FEATURES INTO A MATRIX
%Healthy_Float subject: matrix(9*29)

all_sensors={'LANK', 'LTOE', 'LHIP','LKNE','RANK', 'RTOE', 'RHIP', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

Features.CadenceL=[(120./(diff(PosIC_left.T_01)/100)),(120./(diff(PosIC_left.T_02)/100)),(120./(diff(PosIC_left.T_03)/100))]';
Features.CadenceR=[120./(diff(PosIC_right.T_01)/100),120./(diff(PosIC_right.T_02)/100),120./(diff(PosIC_right.T_03)/100)]';

l=0;
for k=1:length(trials)
    for j=1:length(PosIC_left.(trials{k}))-1
            
   
Features.StanceL(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.Stance.(all_sensors{2});
Features.StanceR(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.Stance.(all_sensors{6});

Features.SwingL(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.Swing.(all_sensors{2});
Features.SwingR(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.Swing.(all_sensors{6});

Features.StepHeightL(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.StepHeight.(all_sensors{2});
Features.StepHeightR(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.StepHeight.(all_sensors{6});

Features.KneeHeightL(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.KneeHeight.(all_sensors{4});
Features.KneeHeightR(j+l,1)=S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.KneeHeight.(all_sensors{8});

Features.HipAngleLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.HipAngle);
Features.HipAngleRmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.HipAngle);

Features.HipAngleLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.HipAngle);
Features.HipAngleRmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.HipAngle);

Features.JointAngleKneeLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.KneeJointAngle);
Features.JointAngleKneeRmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.KneeJointAngle);

Features.JointAngleKneeLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.KneeJointAngle);
Features.JointAngleKneeRmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.KneeJointAngle);

Features.ElevationKneeLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.KneeElevAngle);
Features.ElevationKneeRmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.KneeElevAngle);

Features.ElevationKneeLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.KneeElevAngle);
Features.ElevationKneeLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.KneeElevAngle);

Features.JointAngleAnkleLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.AnkleJointAngle);
Features.JointAngleAnkleRmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.AnkleJointAngle);

Features.JointAngleAnkleLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.AnkleJointAngle);
Features.JointAngleAnkleRmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.AnkleJointAngle);

Features.JointAngleAnkleLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.AnkleElevAngle);
Features.JointAngleAnkleLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.AnkleElevAngle);

Features.JointAngleAnkleLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.AnkleElevAngle);
Features.JointAngleAnkleLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.AnkleElevAngle);

Features.ElevationToeLmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.ToeElevAngle);
Features.ElevationToeRmax(j+l,1)=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.ToeElevAngle);

Features.ElevationToeLmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.LeftAngles.ToeElevAngle);
Features.ElevationToeRmin(j+l,1)=min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).Kin.RightAngles.ToeElevAngle);


    end
 l=l+j;   
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
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsLeft{i})=S6_FLOAT.(trials{k}).(EMGType{l}).(EMGSensorsLeft{i})((10*PosIC_left.(trials{k})(j)):(10*PosIC_left.(trials{k})(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(PosIC_left.(trials{k}))-1
         for i=1:length(EMGSensorsLeft)
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsLeft{i})=S6_FLOAT.(trials{k}).Raw.EMG.(EMGSensorsLeft{i})((10*PosIC_left.(trials{k})(j)):(10*PosIC_left.(trials{k})(j+1)-1));
         end
     end
end




%divide into gait cycles - right
for k=1:length(trials)
     for j=1:length(PosIC_right.(trials{k}))-1
        for l=1:length(EMGType)   
            for i=1:length(EMGSensorsRight)
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.(EMGType{l}).(EMGSensorsRight{i})=S6_FLOAT.(trials{k}).(EMGType{l}).(EMGSensorsRight{i})((10*PosIC_right.(trials{k})(j)):(10*PosIC_right.(trials{k})(j+1)-1));
            end                
         end
     end
end
for k=1:length(trials)
     for j=1:length(PosIC_right.(trials{k}))-1
         for i=1:length(EMGSensorsRight)
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensorsRight{i})=S6_FLOAT.(trials{k}).Raw.EMG.(EMGSensorsRight{i})((10*PosIC_right.(trials{k})(j)):(10*PosIC_right.(trials{k})(j+1)-1));
         end
     end
end


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_FLOAT.(trials{k}).Raw.EMG.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Raw signal');


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_FLOAT.(trials{k}).Filtered.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Filtered signal');


figure
for k=1:length(trials)
     for i=1:4
         subplot(3,4,4*(k-1)+i)
         plot(S6_FLOAT.(trials{k}).Filtered4.(EMGSensors{i}));
         title(sprintf('Trial %d - Sensor %s',k,EMGSensors{i}));
     end
end
suptitle('Final filtered signal');


%--
figure
     for j=2:3
         for i=1:length(fieldnames(S6_FLOAT.(trials{1}).GaitCycles))
            subplot(4,4,4*(j-1)+i)
            plot(S6_FLOAT.(trials{1}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
         end   
     end
suptitle('Trial 1 - Raw');


figure
     for j=2:3
         for i=1:length(fieldnames(S6_FLOAT.(trials{3}).GaitCycles))
            subplot(4,4,4*(j-1)+i)
            plot(S6_FLOAT.(trials{3}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
            title(sprintf('Gait cycle %d - Sensor %s',i,EMGSensors{j}));
         end   
     end
suptitle('Trial 3 - Final filtered');

%% BURSTS CALCULATION FROM STD
%---LTA
for k=1:length(trials)
    for i=1:length(fieldnames(S6_FLOAT.(trials{k}).GaitCycles))
         for j=2%[2,3]
            %Assumption#1.1: for LTA between 750 and 1250 it is always noise
            %Assumption#1.2: we put a higher limit for 2*std
            clear StdNoise;
            clear indeces;
            StdNoise=min([std(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j})(750:1250)),0.0375/2]);
            %StdNoise=max([StdNoise,0.01/2]);
            indeces=find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j})>=(2*StdNoise));
            %Assumption#2: minimum burst length is 350
            Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=[];
            for l=1:(length(indeces)-1)
               if (indeces(l+1)-indeces(l)>350)  
                  Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=cat(1,Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}),indeces(l),indeces(l+1)); 
               end
            end
            
            %Particular case#1.a: no initial (i.e. let's say before first 25% of gait cycle) burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)>0.25*length((S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=1;
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=indeces(1);
            end    
            %Particular case#1.b: no final burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)<0.30*length((S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                 Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=indeces(length(indeces));
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=length(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
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
    for i=1:length(fieldnames(S6_FLOAT.(trials{k}).GaitCycles))
         for j=3%[2,3]
            %Assumption#1.1: for LTA between 750 and 1250 it is always noise
            %Assumption#1.2: we put a higher limit for 2*std
            clear StdNoise;
            clear indeces;
            StdNoise=min([std(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Raw.(EMGSensors{j})(750:1250)),0.0375/2]);
            indeces=find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j})>=(2*StdNoise));
            %Assumption#2: minimum burst length is 350
            Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=[];
            for l=1:(length(indeces)-1)
               if (indeces(l+1)-indeces(l)>350)  
                  Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})=cat(1,Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j}),indeces(l),indeces(l+1)); 
               end
            end
            
            
            %Particular case#1.a: no initial (i.e. let's say before first 25% of gait cycle) burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)>0.25*length((S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=1;
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=indeces(1);
            end    
            %Particular case#1.b: no final burst
            if (isempty(Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})))&&(indeces(1)<0.30*length((S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}))))
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(1,1)=indeces(length(indeces));
                Bursts.Position.(trials{k}).GaitCycles.(numbers{i}).(EMGSensors{j})(2,1)=length(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).EMG.Filtered4.(EMGSensors{j}));
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
%%Bursts calculation visually
% Bursts.Position.T_01.GaitCycles.One.LTA=[203 1353];
% Bursts.Position.T_01.GaitCycles.Two.LTA=[326 1793];
% Bursts.Position.T_01.GaitCycles.Three.LTA=[420 1621];
% Bursts.Position.T_01.GaitCycles.Four.LTA=[250 1588];
% 
% Bursts.Position.T_01.GaitCycles.One.RMG=[563 2038];
% Bursts.Position.T_01.GaitCycles.Two.RMG=[560 1956];
% Bursts.Position.T_01.GaitCycles.Three.RMG=[390 1708];
% Bursts.Position.T_01.GaitCycles.Four.RMG=[410 1773];
% 
% Bursts.Position.T_02.GaitCycles.One.LTA=[333 1650];
% Bursts.Position.T_02.GaitCycles.Two.LTA=[380 1837];
% Bursts.Position.T_02.GaitCycles.Three.LTA=[292 1429];
% 
% Bursts.Position.T_02.GaitCycles.One.RMG=[337 1669];
% Bursts.Position.T_02.GaitCycles.Two.RMG=[453 1700];
% Bursts.Position.T_02.GaitCycles.Three.RMG=[500 1791];
% 
% 
% Bursts.Position.T_03.GaitCycles.One.LTA=[390 1805];
% Bursts.Position.T_03.GaitCycles.Two.LTA=[260 1773];
% 
% Bursts.Position.T_03.GaitCycles.One.RMG=[535 1900];
% Bursts.Position.T_03.GaitCycles.Two.RMG=[458 1805];
% 
%%Calculation of EMG parameters
for i=1:3
    for j=1:length(fieldnames(S6_FLOAT.(trials{i}).GaitCycles))
        for k=2:3
            %Each burst falls between two gait cycles, i.e. each gait
            %cycles has two half-bursts (one at the beginning and one at
            %the end)
            end1burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(1);%end of the first half-burst
            start2burst=Bursts.Position.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})(2);%start of the second half-burst
            lengthburst=length(S6_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Raw.(EMGSensors{k}));
        
            %Duration in seconds (i.e. divided by 1000), i.e. how much of the gait cycle is burst
            Bursts.Duration.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=0.001*(end1burst-1+lengthburst-start2burst);
        
            %Max value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MaxValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=max(S6_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k}));
            
            %Mean value of the burst (in the final filtered signal) in the gait cycle
            Bursts.MeanValue.(trials{i}).GaitCycles.(numbers{j}).(EMGSensors{k})=mean(S6_FLOAT.(trials{i}).GaitCycles.(numbers{j}).EMG.Filtered4.(EMGSensors{k})([1:end1burst,start2burst:lengthburst]));
        end
    end
end

%% MATRIC FOR FEATURES PCA

l=0;
for k=1:length(trials)
    for j=1:length(PosIC_left.(trials{k}))-1
           

Features.DurationBurstLTA(j+l,1)=Bursts.Duration.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{2});
Features.DurationBurstRGM(j+l,1)=Bursts.Duration.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{3});
Features.MaxAmplitudeBurstLTA(j+l,1)=Bursts.MaxValue.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{2});
Features.MaxAmplitudeBurstRGM(j+l,1)=Bursts.MaxValue.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{3});
Features.MeanAmplitudeBurstLTA(j+l,1)=Bursts.MeanValue.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{2});
Features.MeanAmplitudeBurstRGM(j+l,1)=Bursts.MeanValue.(trials{k}).GaitCycles.(numbers{j}).(EMGSensors{3});
        
    end
    l=l+j;
end