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



%% Divide in gait cycles

