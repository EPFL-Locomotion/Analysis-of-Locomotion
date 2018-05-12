clear all
close all
clc

% adding the paths and loading data

addpath(genpath('Healthy Recordings\Subject6\FLOAT'));

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
    %pause(0.0001);
    
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
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).StepHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3));
            
        end
        
    end
end

% max height of the knee

sensors={'LKNE', 'RKNE'};

for k=1:length(trials)
    
    for i=1:length(numbers)
        
        for j=1:length(sensors)
         
            % max knee height
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).KneeHeight.(sensors{j})=max(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3));
            
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
        
        for j=1:length(numbers)
            
            for i=1:size(S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,1}),1)
                
                TOE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,1})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,1})(i,3)];
                ANK_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,2})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,2})(i,3)];
                KNE_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,3})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,3})(i,3)];
                HIP_points(i,:)=[S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,4})(i,2) S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,4})(i,3)];
                
                VectorANKLE_TOE(i,:)=[(ANK_points(i,1)-TOE_points(i,1)) (ANK_points(i,2)-TOE_points(i,2))];
                VectorKNEE_ANKLE(i,:)=[(KNE_points(i,1)-ANK_points(i,1)) (KNE_points(i,2)-ANK_points(i,2))];
                VectorHIP_KNEE(i,:)=[(HIP_points(i,1)-KNE_points(i,1)) (HIP_points(i,2)-KNE_points(i,2))];
                
                Vertical = [0 1];
                Horizontal = [1 0];
                
                % hip angle to get extension/flexion
                Angle_hip_vertical(i)=acos(dot(VectorHIP_KNEE(i,:),Vertical)/(norm(VectorHIP_KNEE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).HipAngle(i)=Angle_hip_vertical(i)*180/pi;
                
                % joint angle knee
                Angle_knee(i)=acos(dot(VectorHIP_KNEE(i,:),VectorKNEE_ANKLE(i,:))/(norm(VectorHIP_KNEE(i,:))*norm(VectorKNEE_ANKLE(i,:))));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).KneeJointAngle(i)=Angle_knee(i)*180/pi;
                
                % elevation angle of the knee joint
                Angle_knee_vertical(i)=acos(dot(VectorKNEE_ANKLE(i,:),Vertical)/(norm(VectorKNEE_ANKLE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).KneeElevAngle(i)=Angle_knee_vertical(i)*180/pi;
                
                % joint angle ankle
                Angle_ankle(i)=acos(dot(VectorKNEE_ANKLE(i,:),VectorANKLE_TOE(i,:))/(norm(VectorKNEE_ANKLE(i,:))*norm(VectorANKLE_TOE(i,:))));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).AnkleJointAngle(i)=Angle_ankle(i)*180/pi;
                
                % elevation angle af the ankle joint
                Angle_ankle_vertical(i)=acos(dot(VectorANKLE_TOE(i,:),Vertical)/(norm(VectorANKLE_TOE(i,:))*norm(Vertical)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).AnkleElevAngle(i)=Angle_ankle_vertical(i)*180/pi;
                
                % elevation toe from the ground
                Angle_toe_hor(i)=acos(dot(VectorANKLE_TOE(i,:),Horizontal)/(norm(VectorANKLE_TOE(i,:))*norm(Horizontal)));
                S6_FLOAT.(trials{k}).GaitCycles.(numbers{j}).(sensors{n,5}).ToeElevAngle(i)=Angle_toe_hor(i)*180/pi;
                
            end
        end
    end
end




