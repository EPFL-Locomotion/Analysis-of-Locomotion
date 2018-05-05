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

IC_right(1)=min(S6_FLOAT.T_01.Raw.Kin.RANK(70:250,3));
IC_right(2)=min(S6_FLOAT.T_01.Raw.Kin.RANK(251:500,3));
IC_right(3)=min(S6_FLOAT.T_01.Raw.Kin.RANK(501:750,3));
IC_right(4)=min(S6_FLOAT.T_01.Raw.Kin.RANK(751:end,3));

PosIC_right(5)=length(S6_FLOAT.T_01.Raw.Kin.RANK);

for i=1:4
    
PosIC_right(i)=find(S6_FLOAT.T_01.Raw.Kin.RANK(:,3)==IC_right(i));

end

%% split in cycles left

IC_left(1)=min(S6_FLOAT.T_01.Raw.Kin.LANK(70:250,3));
IC_left(2)=min(S6_FLOAT.T_01.Raw.Kin.LANK(251:500,3));
IC_left(3)=min(S6_FLOAT.T_01.Raw.Kin.LANK(501:750,3));
IC_left(4)=min(S6_FLOAT.T_01.Raw.Kin.LANK(751:end,3));

PosIC_left(5)=length(S6_FLOAT.T_01.Raw.Kin.LANK);

for i=1:4
    
PosIC_left(i)=find(S6_FLOAT.T_01.Raw.Kin.LANK(:,3)==IC_left(i));

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
plot(S6_FLOAT.T_01.GaitCycles.One.RANK(:,2), S6_FLOAT.T_01.GaitCycles.One.RANK(:,3)); 

figure;
plot(S6_FLOAT.T_01.GaitCycles.Two.RTOE(:,2), S6_FLOAT.T_01.GaitCycles.Two.RTOE(:,3)); 



%% search of parameters

sensors={ 'LTOE', 'RTOE'};
numbers={'One', 'Two', 'Three', 'Four'};
trials={'T_01', 'T_02', 'T_03'};

% stance duration (divide by fsKin to get it in seconds)

for k=1:length(trials)
    
    for i=1:length(numbers)
        
        for j=1:length(sensors)
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Stance.(sensors{j})=find(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3)==min(S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j})(:,3)))\100; 
        
        end
        
    end
end

% swing duration

sensors={ 'LANK', 'RANK'};

for k=1:length(trials)
    
    for i=1:length(numbers)
        
        for j=1:length(sensors)
            A=S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).(sensors{j});
            S6_FLOAT.(trials{k}).GaitCycles.(numbers{i}).Swing.(sensors{j})=(find(A(:,3)==max(A(:,3)))-find(A(:,3)==min(A(1:61,3))))/100;

        end
        
    end
end



% step height

% angles




