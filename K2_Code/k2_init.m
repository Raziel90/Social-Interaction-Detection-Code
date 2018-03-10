clear
%basepath = '/media/claudio/52FFC5D351F60A1C/';
%basepath='/media/ccoppola/52FFC5D351F60A1C/';
%'/media/ccoppola/52FFC5D351F60A1C'

folpath = fopen('config');
folders = textscan(folpath,'%s');
fclose(folpath);
codefolder = folders{1}{1};
datafolder = folders{1}{2};
annotationpath = folders{1}{3};
outputfolder = folders{1}{4};
status = rmdir([outputfolder,'/*']);
mkdir(outputfolder,'k2_figures');

clear folders folpath
%depthpath=[basepath,'/Social_Activity_Detection_Data/Output/'];
%addpath(genpath([basepath,'/Social-Activities-Code']));
%addpath([basepath,'/Social_Activity_Detection_Data']);
depthpath = datafolder;
addpath(genpath(codefolder))
addpath(genpath(datafolder))
%clear basepath;
% load('dataset.mat')
k2_Import_videos_annotations
% [Glo,Det,Denerg,DistOrient]=k2_segmentation_feat(skeleton_bags,[],basepath,0);