function [qtc_sequences,qtc_distances] = k2_convert_qtc_all(skeleton_bags,QTC_threshold,joint_idx,convert2indexflag)

qtc_sequences=cell(size(skeleton_bags,1),1);
qtc_distances=cell(size(skeleton_bags,1),1);
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
    [q1,q2,q3,q4,qd]=QTCC(squeeze(skel1(1,1,:))',squeeze(skel1(2,1,:))',squeeze(skel2(1,1,:))',squeeze(skel2(2,1,:))');
    qtccnum=[q1;q2;q3;q4]';
%     [qtccnum]=qtcc(squeeze(skel1(:,:))',squeeze(skel2(:,:))',QTC_threshold,1,0);
    
    if convert2indexflag ==1 
        qtc_sequences{k}=qtc2idx([qtccnum([1,1],:);qtccnum]);
    else
        qtc_sequences{k}=[qtccnum([1,1],:);qtccnum];
    end
    qtc_distances{k}=[0;0;qd'];
end


end