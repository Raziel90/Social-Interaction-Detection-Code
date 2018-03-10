function [ output_args ] = k2_load_depth( basepath )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%bp='/media/ccoppola/52FFC5D351F60A1C/Social_Activity_Detection_Data/Output/'
recordingnames=dir([basepath,'/session_*/recording_*']);
DepthVideos=repmat(struct('session','','bagname','','data',[]),length(recordingnames),1);
for k=1:length(recordingnames)
    DepthVideos(k).bagname=string(recordingnames(k).name);
    splitted=strsplit(recordingnames(k).folder);
    DepthVideos(k).session=string(splitted(end));
   depthfiles=dir([recordingnames(k).folder,'/',recordingnames(k).name,'/depth/*.txt']); 
   filepath=strcat({depthfiles.folder},'/',{depthfiles.name});
   
   depthvideo=cell(length(filepath),1);
   for n=1:length(filepath)
       depthvideo{n}=load(filepath{n});
   end
   DepthVideos(k).data=cell2mat(depthvideo);
   
   
end

end

