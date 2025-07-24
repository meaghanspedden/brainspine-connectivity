%% main script
clear all
close all

%% give sub as input here?
%also need to automatise that file name for 212rest is 002

%% after running preproc
tic
sub='OP00215';
analysis='static';

addpath('D:\brainspineconnectivity')
setup_source

create_leadfields

beamforming

beamforming_brain

source_freq

toc
