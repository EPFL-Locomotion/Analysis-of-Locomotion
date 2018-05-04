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

load('SCI_NOFLOAT_T02.m');
load('SCI_NOFLOAT_T01.m');
load('S6_NOFLOAT_T01.m');
load('NO_FLOAT_CRUTCHES.mat');
load('FLOAT_NO_CRUTCHES.mat');
load('S6_NO_FLOAT.mat');
load('S6_FLOAT.mat');

%% divide in gait cycle

%% emg analysis filtering the data