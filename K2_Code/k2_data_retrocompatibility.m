function [ Sessions ] = k2_data_retrocompatibility( skeleton_bags )
%k2_data_retrocompatibility converts the data standard of the Continous RGB-D
%Social Activity Dataset to the one used in Social Activity Dataset (IROS2016)
%   Detailed explanation goes here

    %conversion table between new and old skeleton
    converted_joints=[5,4,2,6,7,10,11,14,15,18,19,9,12,16,20];
    
    newnumdimensions=7;
    newnumjoints=25;
    oldnumdimensions=6;
    numbags=size(skeleton_bags,1);

    Sessions=repmat(struct('action',{[]}),numbags,2);



    for k=1:numbags
        datalen1=length(skeleton_bags{k,3}(1).data);
        datalen2=length(skeleton_bags{k,3}(2).data);
        % lenclip=length(skeleton_bags{k,2});

        converted_k2_reshaped_skel1=zeros(oldnumdimensions,length(converted_joints),datalen1);
        converted_k2_reshaped_skel2=zeros(oldnumdimensions,length(converted_joints),datalen2);
        k2_reshaped_skel1=nan(newnumdimensions,newnumjoints,datalen1);
        k2_reshaped_skel2=nan(newnumdimensions,newnumjoints,datalen2);
        k2_reshaped_skel1(:,:,1:datalen1)=reshape(skeleton_bags{k,3}(1).data',newnumdimensions,newnumjoints,datalen1);
        k2_reshaped_skel2(:,:,1:datalen2)=reshape(skeleton_bags{k,3}(2).data',newnumdimensions,newnumjoints,datalen2);
        converted_k2_reshaped_skel1(1:3,:,:)=k2_reshaped_skel1([1,3,2],converted_joints,:);
        converted_k2_reshaped_skel2(1:3,:,:)=k2_reshaped_skel2([1,3,2],converted_joints,:);

        for j=1:min(datalen1,datalen2)
            [converted_k2_reshaped_skel1([4],:,j),converted_k2_reshaped_skel1([5],:,j),converted_k2_reshaped_skel1([6],:,j)]=quat2angle(k2_reshaped_skel1([7,4,5,6],converted_joints,j)', 'XYZ');
            [converted_k2_reshaped_skel2([4],:,j),converted_k2_reshaped_skel2([5],:,j),converted_k2_reshaped_skel2([6],:,j)]=quat2angle(k2_reshaped_skel1([7,4,5,6],converted_joints,j)', 'XYZ');
        end
        Sessions(k,1).action{1}=reshape(converted_k2_reshaped_skel1,oldnumdimensions*length(converted_joints),datalen1)';
        Sessions(k,2).action{1}=reshape(converted_k2_reshaped_skel2,oldnumdimensions*length(converted_joints),datalen2)';
%         skeleton_bags{k,3}(2).temporal_presence
        mutual_presence=skeleton_bags{k,3}(1).temporal_presence&skeleton_bags{k,3}(2).temporal_presence;
        Sessions(k,2).action{1}(min(datalen1,datalen2)+1:end,:)=[];
        Sessions(k,1).action{1}(min(datalen1,datalen2)+1:end,:)=[];
%         Sessions(k,1).action{1}=Sessions(k,1).action{1}(mutual_presence,:);
%         Sessions(k,1).action{1}(~mutual_presence(1:datalen1),:)=[];
%         Sessions(k,2).action{1}=Sessions(k,2).action{1}(mutual_presence,:);
%         Sessions(k,2).action{1}(~mutual_presence(1:datalen2),:)=[];
    end

end

