function [ skeleton_bags ] = k2_skeleton_text_to_matlab( base_path )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
base_path='/home/ccoppola/Social_Activity_Detection_Data/Output';
end
numfol=length(dir([base_path,'/session_*']));
for h=1:numfol
skeleton_files(h,:)=dir([base_path,'/session*_',num2str(h),'/**/skeletons.txt']);
end
skeleton_files=reshape(skeleton_files',numel(skeleton_files),1);
skeleton_bags=cell(length(skeleton_files),3);
for i=1:size(skeleton_files,1)
%     skelfile=fopen('/home/ccoppola/Social_Activity_Detection_Data/Output/session_1/recording_1482159915/skeletons.txt','r');
    skelfile=fopen([skeleton_files(i).folder,'/',skeleton_files(i).name],'r');
    timed_Data=textscan(skelfile,'%f,%s\n','Delimiter','');
%     time=vertcat(timed_Data{:,1}-min(timed_Data{:,1}));
    time=vertcat(timed_Data{:,1}-min(timed_Data{:,1}));
    fclose(skelfile);
    string_Data=split(string(timed_Data{2}),',');
    numerical_rows=cellfun(@(X)str2num(X),cellstr(string_Data),'UniformOutput',false);
    empty_rows=cell2mat(cellfun(@(X)isempty(X),numerical_rows,'UniformOutput',false));
    nonempty_cols=sum(empty_rows)~=size(numerical_rows,1);
    user_skeletons=repmat(struct('trackingId','','temporal_presence',false(size(empty_rows,1),1),'data',zeros), sum(nonempty_cols), 1 );
    for skel_num=1:length(user_skeletons)
        mat=vertcat(numerical_rows{~empty_rows(:,skel_num),skel_num});
        user_skeletons(skel_num).trackingId=cell2mat(unique(cellfun(@strtok,cellstr(string_Data(~empty_rows(:,skel_num),skel_num)),'UniformOutput',false))); %from the numerical part i lost it, so i take it from the string
        user_skeletons(skel_num).temporal_presence(~empty_rows(:,skel_num))=true;
        user_skeletons(skel_num).data=mat(:,2:end);
    end
    splitted= strsplit('/',skeleton_files(i).folder);
    skeleton_bags{i,1}=string(splitted(end));
    skeleton_bags{i,2}=time;
    skeleton_bags{i,3}=user_skeletons;
end


end

