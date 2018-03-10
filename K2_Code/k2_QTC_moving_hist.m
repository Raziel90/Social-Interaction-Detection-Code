function [ qtc_hists ] = k2_QTC_moving_hist(skeleton_bags,window,QTC_threshold,joint_idx,doublehistflag)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




qtc_sequences=cell(size(skeleton_bags,1),1);
qtc_hists=cell(size(skeleton_bags,1),1);
qtcdistk=cell(size(skeleton_bags,1),1);
numdimensions=7;
numjoints=25;
for k=1:size(skeleton_bags,1)
    
    lenclip=length(skeleton_bags{k,2});
    skel1=nan(numdimensions,numjoints,lenclip);
    skel2=nan(numdimensions,numjoints,lenclip);
    skel1(:,:,skeleton_bags{k,3}(1).temporal_presence)=reshape(skeleton_bags{k,3}(1).data',numdimensions,numjoints,length(skeleton_bags{k,3}(1).data));
    skel2(:,:,skeleton_bags{k,3}(2).temporal_presence)=reshape(skeleton_bags{k,3}(2).data',numdimensions,numjoints,length(skeleton_bags{k,3}(2).data));
    skel1=skel1([1,3],joint_idx,:);
    skel2=skel2([1,3],joint_idx,:);
%     [qtccnum]=qtcc(squeeze(skel1(:,:))',squeeze(skel2(:,:))',QTC_threshold,1,0);
    [q1,q2,q3,q4,qd]=QTCC(squeeze(skel1(1,1,:))',squeeze(skel1(2,1,:))',squeeze(skel2(1,1,:))',squeeze(skel2(2,1,:))');
    qtccnum=[q1;q2;q3;q4]';
    if doublehistflag ==1 
        qtc_sequences{k}=[qtc2idx(qtccnum(:,1:2),2),qtc2idx(qtccnum(:,3:4),2)];
        qtc_hists{k}=[hist(im2col([zeros(window+1,1);qtc_sequences{k}(:,1)],[window,1],'sliding'),9),hist(im2col([zeros(window-1,1);qtc_sequences{k}(:,1)],[window,1],'sliding'),9)]';
        
    else
        qtc_sequences{k}=qtc2idx(qtccnum,4);
        qtc_hists{k}=hist(im2col([zeros(window+1,1);qtc_sequences{k}(:,1)],[window,1],'sliding'),81)';
    end
    qtcdistk{k}=qd';
    qtc_hists{k}(:,sum(qtc_hists{k})~=0)=qtc_hists{k}(:,sum(qtc_hists{k})~=0)./sum(qtc_hists{k}(:,sum(qtc_hists{k})~=0));
end


end

