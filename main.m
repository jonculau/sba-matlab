close all;
clear all;
clc;

%@Author: Jonathan Culau
%% id
r = 1; % rotational id
t = 2; % translational id
g = 3; % unified
addpath('objects');
addpath('utils');

%% Initialization
sp = StewartPlatform;

ControlDesign;
Simulation;

