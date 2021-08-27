%% This file is used to create the geometry (hex mesh) for the ML dataset
% Pipeline: Tree skeletion -> Tree hex mesh -> Extract pipe/bif mesh

%% (Must run) Skeleton extraction for training and prediction
% Create the skeleton information of each pipe and bifurcation
% The skeleton information is then used for extract the mesh of pipe or bif
clc;
clear;
addpath(genpath(pwd));
start_trees;

%% Generate quad mesh
clc;
clear;
addpath(genpath(pwd));
start_trees;

io_path = '..//io//single_pipe//';

GenHex_ExtractQuad(io_path);
