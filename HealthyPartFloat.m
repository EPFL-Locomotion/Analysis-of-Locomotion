clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Separating into gait cycles'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\NO_FLOAT_CRUTCHES\GAIT FILES'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\MAT'));
addpath(genpath('SCI subject (updated file)\FLOAT_NO_CRUTCHES\GAIT FILES'));
addpath(genpath('Healthy Recordings\Subject6\NO_FLOAT'));
addpath(genpath('Healthy Recordings\Subject6\FLOAT'));
addpath(fullfile('Separating into gait cycles')); 

load('NO_FLOAT_CRUTCHES.mat');
load('FLOAT_NO_CRUTCHES.mat');
load('S6_NO_FLOAT.mat');
load('S6_FLOAT.mat');


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
    pause(0.01);
    
end
hold off;

%% split in cycles right 

IC_right(1)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(1:250,3));
IC_right(2)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(251:500,3));
IC_right(3)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(501:700,3));
IC_right(4)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(701:900,3));
IC_right(5)=max(S6_FLOAT.T_01.Raw.Kin.RTOE(901:end,3));

for i=1:5
    
PosIC_right(i)=find(S6_FLOAT.T_01.Raw.Kin.RTOE(:,3)==IC_right(i));

end

%% split in cycles left

IC_left(1)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(1:50,3));
IC_left(2)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(51:400,3));
IC_left(3)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(401:600,3));
IC_left(4)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(601:800,3));
IC_left(5)=max(S6_FLOAT.T_01.Raw.Kin.LTOE(800:end,3));


for i=1:5
    
PosIC_left(i)=find(S6_FLOAT.T_01.Raw.Kin.LTOE(:,3)==IC_left(i));

end

%% making struct for each cycles and sensor

sensors_left={'LANK', 'LTOE', 'LHIP','LKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for i=1:length(sensors_left)
        for j=1:length(numbers)
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors_left{i})=S6_FLOAT.T_01.Raw.Kin.(sensors_left{i})(PosIC_left(j):PosIC_left(j+1),:);
        end
    end
end

sensors_right={'RANK', 'RTOE', 'RHIP', 'RKNE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    for i=1:length(sensors_right)
        for j=1:length(numbers)
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors_right{i})=S6_FLOAT.T_01.Raw.Kin.(sensors_right{i})(PosIC_right(j):PosIC_right(j+1),:);
        end
    end
end

%% sample plots 
%(avoiding to plot the x direction cause we suppose he moves straigh)

figure;
plot(S6_FLOAT.T_01.GaitCycles.One.LANK(:,2), S6_FLOAT.T_01.GaitCycles.One.LANK(:,3)); 

figure;
plot(S6_FLOAT.T_01.GaitCycles.Two.RTOE(:,2), S6_FLOAT.T_01.GaitCycles.Two.RTOE(:,3)); 



%% Distance parameters

sensors={'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

% divide by fsKin to get it in seconds

for k=1:length(trials)
    
    for i=1:length(numbers)
        
        for j=1:length(sensors)
            
            % stance
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Stance.(sensors{j})=find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3)==min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3)))\100; 
            
            % swing
            A=S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j});
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Swing.(sensors{j})=(find(A(:,3)==max(A(150:end,3)))-find(A(:,3)==min(A(100:end,3))))/100;

            % max toe step height
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).StepHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3))
            
        end
        
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(numbers)
        
        for j=1:length(sensors)
         
            % max knee height
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).KneeHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3))
            
        end
    end
end



%% angles for LEFT foot

numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for j=1:length(numbers)
        
        for i=1:size(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LTOE,1)
            
            LTOE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LTOE(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LTOE(i,3)];
            LANK_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LANK(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LANK(i,3)];
            LKNE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LKNE(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LKNE(i,3)];
            LHIP_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LHIP(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LHIP(i,3)];
            
            VectorANKLE_TOE(i,:)=[(LANK_points(i,1)-LTOE_points(i,1)) (LANK_points(i,2)-LTOE_points(i,2))];
            VectorKNEE_ANKLE(i,:)=[(LKNE_points(i,1)-LANK_points(i,1)) (LKNE_points(i,2)-LANK_points(i,2))];
            VectorHIP_KNEE(i,:)=[(LHIP_points(i,1)-LKNE_points(i,1)) (LHIP_points(i,2)-LKNE_points(i,2))];
            
            Vertical = [0 1];
            Horizontal = [1 0];
            
            % hip angle to get extension/flexion
            Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.HipAngle(i)=Angle_hip_vertical(i)*180/pi;
            
            % joint angle knee
            Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.KneeJointAngle(i)=Angle_knee(i)*180/pi;
            
            % elevation angle of the knee joint
            Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
            
            % joint angle ankle
            Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
            
            % elevation angle af the ankle joint
            Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
            
            % elevation toe from the ground
            Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).LeftAngles.ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
            
        end
    end
end

%% angles for RIGHT foot

numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

for k=1:length(trials)
    
    for j=1:length(numbers)
        
        for i=1:size(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RTOE,1)
            
            RTOE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RTOE(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RTOE(i,3)];
            RANK_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RANK(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RANK(i,3)];
            RKNE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RKNE(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RKNE(i,3)];
            RHIP_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RHIP(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RHIP(i,3)];
            
            VectorANKLE_TOE(i,:)=[(RANK_points(i,1)-RTOE_points(i,1)) (RANK_points(i,2)-RTOE_points(i,2))];
            VectorKNEE_ANKLE(i,:)=[(RKNE_points(i,1)-RANK_points(i,1)) (RKNE_points(i,2)-RANK_points(i,2))];
            VectorHIP_KNEE(i,:)=[(RHIP_points(i,1)-RKNE_points(i,1)) (RHIP_points(i,2)-RKNE_points(i,2))];
            
            Vertical = [0 1];
            Horizontal = [1 0];
            
            % angle to get extension/flexion
            Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.HipAngle(i)=Angle_hip_vertical(i)*180/pi;
            
            % joint angle knee
            Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.KneeJointAngle(i)=Angle_knee(i)*180/pi;
            
            % elevation angle of the knee joint
            Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
            
            % joint angle ankle
            Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
            
            % elevation angle af the ankle joint
            Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
            
            % elevation toe from the ground
            Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).RightAngles.ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
            
        end
    end
end



