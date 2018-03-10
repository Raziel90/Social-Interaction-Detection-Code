function [... %{depenerg, movingpixels,variation,%} 
    gloentr,detentr,detenerg,uorient... %{,cumulative_joint_energ
    ] = segmentation_feat(Sessions,Depthclips,basepath,plotit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


numplots=8;
fps=15;
w=1;
th1=0.2; 
dispth=0.2;
order=2;        
    
    


gloentr=vel_entropy(Sessions,fps,w);

detentr=vel_entropy_detailed(Sessions,fps,w);
detenerg=energy_detailed(Sessions,fps,w);
[uorient,Fovdist,Fpdist]=Ubodyorientdist(Sessions,fps,w);
% [depenerg,movingpixels,variation]=Depthenerg(Depthclips,15,w,4,dispth);
cumulative= minmax_normalization(accumulate_joint_energy(detenerg,w,[1:2,4:15]));
Joint_str={'SpineBase', 'SpineMid', 'SpineShoulder', 'Neck', 'Head',...
'ShoulderLeft', 'ElbowLeft', 'WristLeft', 'HandLeft','ShoulderRight',...
'ElbowRight', 'WristRight', 'HandRight','HipLeft', 'KneeLeft', 'AnkleLeft',...
'FootLeft','HipRight', 'KneeRight', 'AnkleRight', 'FootRight','HandTipLeft', ...
'ThumbLeft','HandTipRight', 'ThumbRight'};



outsignal=cell(1,size(Sessions,1));
thresholded=cell(1,size(Sessions,1));

%figure(1),suptitle({['normalized cumulative energy '],['']})
%figure(2),suptitle({['thresholded th= ',num2str(th1)],['']})
%figure(3),suptitle({['th= ',num2str(th1),' filtered order= ',num2str(order)],['']})
for c=1:size(Sessions,1)
    
    
    
    
    thresholded{c}=cumulative{c}>th1;
    
    if plotit==1
        
        T=fft(thresholded{c});
        [~,sorted]=sort(T,'descend');
        S=zeros(size(T));
        S(sorted(1:2*order+1))=T(sorted(1:2*order+1));
        outsignal{c}=ifft(S);
        h=figure(c);suptitle({['clip no: ',num2str(c)],['']})
        set(h, 'units', 'centimeters');
        set(h, 'position', [3 3 53 53]);
        
%         subplot(ceil(numplots/4),4,1),plot(w+1:length(depenerg{c}),depenerg{c}(6:end)),hold on,plot(w+1:length(depenerg{c}),std_signal(depenerg{c}(6:end),2*w)), title(['energy computed from depth'])
%         subplot(ceil(numplots/4),4,2),plot(w+1:length(movingpixels{c}),movingpixels{c}(6:end)),hold on,plot(w+1:length(movingpixels{c}),std_signal(movingpixels{c}(6:end),2*w)),title(['number of pixel resulting from the difference'])
        subplot(ceil(numplots/4),4,3),plot(w+1:size(cumulative{c},2)+w,cumulative{c}),hold on,plot(w+1:length(cumulative{c})+w,std_signal(cumulative{c},2*w)),title({['normalized cumulative energy'],[' of joints related to torso']})
        subplot(ceil(numplots/4),4,4),plot(1:size(detenerg{c,1}(3,w+1:end),2),detenerg{c,1}(3,w+1:end)),hold on,plot(w+1:length(detenerg{c}),std_signal(detenerg{c,1}(3,w+1:end),2*w)), title(['absolute energy of torso\_user1'])
        subplot(ceil(numplots/4),4,5),plot(1:size(detenerg{c,2}(3,w+1:end),2),detenerg{c,2}(3,w+1:end)),hold on,plot(w+1:length(detenerg{c}),std_signal(detenerg{c,2}(3,w+1:end),2*w)), title(['absolute energy of torso\_user2'])
        %subplot(ceil(numplots/4),4,6),plot(thresholded{c}),ylim([0 1.1]),title(['thresholded th= ',num2str(th1)])
        %subplot(ceil(numplots/4),4,7),plot(outsignal{c}),ylim([0 1]),title(['filtered order= ',num2str(order)])
        subplot(ceil(numplots/4),4,6),plot(rad2deg(uorient{c}(1,:))),ylim([0 361]),hold on,plot(1:length(uorient{c}(1,:)),std_signal(uorient{c}(1,:),2*w)),title(['angle between orientation of the users'])
        %subplot(ceil(numplots/4),4,7),plot(rad2deg(uorient{c}(3,:)-(uorient{c}(2,:)))),ylim([-360 361]),title(['orientation user\_2'])
        subplot(ceil(numplots/4),4,7),plot(uorient{c}(2,:)),hold on,plot(1:length(uorient{c}(1,:)),std_signal(uorient{c}(2,:),2*w)),title(['distance between torsos'])
%         subplot(ceil(numplots/4),4,8),hold on,plot(Sessions{c,2},double(Sessions{c,3}(1).temporal_presence)),plot(Sessions{c,2},double(Sessions{c,3}(2).temporal_presence)+1),plot(Sessions{c,2},2+double(Sessions{c,3}(1).temporal_presence&Sessions{c,3}(2).temporal_presence)+1),ylim([-0.1 3.1])

%         savefig([basepath,'Social-Activities-Code/k2_figures/','clip-',num2str(c)])
%         saveas(h,[basepath,'Social-Activities-Code/k2_figures/','clip-',num2str(c)],'epsc')
        saveas(h,[basepath,'Social-Activities-Code/k2_figures/','clip-',num2str(c)],'png')
    end
end
% figure(1)
% for i=1:15
% subplot(4,4,i),plot((1:length(detentr{1,1}))/15,detentr{1,1}(i,:)'),title(Joint_str{i}),ylabel('log energy'),ylim([-50 1])
% end
% subplot(4,4,i+1),plot((1:length(detentr{1,1}))/15,gloentr{1,1}),title('global user 1'),ylabel('log energy'),ylim([-500 1])
% 
% 
% 
% for i=1:15
%     df=fps/length(detentr{1,1}(i,6:end));
%     f=-fps/2+df:df:fps/2;
%     figure(2),subplot(4,4,i),plot(f,abs(fftshift(fft(detentr{1,1}(i,6:end)')))),title(Joint_str{i}),ylabel('Amp log energy'),ylim([0 200])
% end
% figure(2),subplot(4,4,i+1),plot(f,abs(fftshift(fft(gloentr{1,1}(6:end))))),title('global user 1'),ylabel('Amp log energy'),ylim([0 1000])
% 
% 
% for i=1:15
% figure(3),subplot(4,4,i),plot(f,angle(fftshift(fft(detentr{1,1}(i,6:end)')))),title(Joint_str{i}),ylabel('Angle log energy')%,ylim([ 1])
% end
% figure(3),subplot(4,4,i+1),plot(f,angle(fftshift(fft(gloentr{1,1}(6:end))))),title('global user 1'),ylabel('Angle log energy')%,ylim([-500 1])
% 


end


function [ Uorientdist,Fovdist,Fpdist ] = Ubodyorientdist(Sessions,fps,w)


    slack=0.3;
    cond_th=0.1;
    
    numjoints=25;
    numdimensions=7;
    numsubjects=2;
%     pos_columns=(reshape([1:7:175;(1:7:174)+1;(1:7:173)+2],1,length([1:7:175,(1:7:174)+1,(1:7:173)+2])));

    %noise_model=@(X)normpdf(X,0,0.01)*0.001;
    
    [numclips]=size(Sessions,1);
    Uorientdist =cell(numclips,1);
    Fovdist =cell(numclips,2);
    Fpdist =cell(numclips,2);
    
    for i=1:numclips
       lenclip=length(Sessions{i,2});
       skel1=nan(numdimensions,numjoints,lenclip);
       skel2=nan(numdimensions,numjoints,lenclip);
       skel1(:,:,Sessions{i,3}(1).temporal_presence)=reshape(Sessions{i,3}(1).data',numdimensions,numjoints,length(Sessions{i,3}(1).data));
       skel2(:,:,Sessions{i,3}(2).temporal_presence)=reshape(Sessions{i,3}(2).data',numdimensions,numjoints,length(Sessions{i,3}(2).data));
       skel1=skel1(1:3,:,:);
       skel2=skel2(1:3,:,:);
       istants_of_common_presence=find(Sessions{i,3}(1).temporal_presence&Sessions{i,3}(2).temporal_presence);
       
       
            tg1=squeeze(skel1([1,2],6,:)-skel1([1,2],10,:));
            tg2=squeeze(skel2([1,2],6,:)-skel2([1,2],10,:));
            n1=[-tg1(2,:);tg1(1,:)];
            n2=[-tg2(2,:);tg2(1,:)];
            d1=dot(n1,squeeze(skel1([1,2],2,:)));
            d2=dot(n2,squeeze(skel2([1,2],2,:)));
            c1=dot(tg1,squeeze(skel1([1,2],2,:)));
            c2=dot(tg2,squeeze(skel2([1,2],2,:)));
            X=zeros(size(n1));
            R=zeros(size(n1,2),1);
            %find interceptions
            for k=1:size(n1,2)
                [X(:,k),R(k)]=linsolve([tg1(:,k)';tg2(:,k)'],[c1(k);c2(k)]);
            end
            Fp_segmentation=(dot(X,n1)-d1>-slack&dot(X,n2)-d2>-slack);
            parallel_idx=R<cond_th;
            Parallel_segmentation=dot(squeeze(skel2([1,3],3,:)),n1)-d1>-slack&dot(squeeze(skel1([1,3],3,:)),n2)-d2>-slack;
            
            %n1_norm=n1/norm(n1);
            %n2_norm=n1/norm(n2);
            %angle=abs(atand(n1(2,:)./n1(1,:))-atand(n2(2,:)./n1(1,:)));
            %angle=mod(atan2(n1(2,:),n1(1,:))-atan2(n2(2,:),n2(1,:)),pi);
            cosangle=abs(cos(abs(atan2(n1(2,:),n1(1,:))-atan2(n2(2,:),n2(1,:)))));

%        feat=nan(6,size(skel1,3));
            feat=nan(6,lenclip);
%        for t=1:lenclip
% 
%             %line(x,y,'color','k','LineWidth',2)
%             %normal = [skel1([1],3,t),skel1([3],3,t)] + null(skel1([1,3],4,t)'-skel1([1,3],6,t)')';
%             %normal_mid1 = [mean(skel1(1,[4,6],t)),mean(skel1([3],[4,6],t))] + null(skel1([1,3],4,t)'-skel1([1,3],6,t)')';
%             %normal_mid2 = [mean(skel2(1,[4,6],t)),mean(skel2([3],[4,6],t))] + null(skel2([1,3],4,t)'-skel2([1,3],6,t)')';
%       
%             
% 
%             %line([mean(x),normal(1)],[mean(y),normal(2)],'color','r','LineWidth',2)
%             %sk1l=skel1([1,3],4,t);
%             %sk1r=skel1([1,3],6,t);
%             
% %             feat(1,t)=atan(mean(skel1(3,[4,6],t))-normal_mid1(2)./mean(skel1(1,[4,6],t))-normal_mid1(1));
% %             feat(2,t)=atan(mean(skel2(3,[4,6],t))-normal_mid2(2)./mean(skel2(1,[4,6],t))-normal_mid2(1));
% %             feat(3,t)=atan2(skel1([3],3,t)-skel2([3],3,t),skel1([1],3,t)-skel2([1],3,t));
% %             feat(4,t)=sqrt(sum((skel1([1,3],3,t)-skel2([1,3],3,t)).^2));
% 
%             feat(1,:)=angle;
%             feat(2,t)=sqrt(sum((skel1([1,3],3,t)-skel2([1,3],3,t)).^2));
%             
% 
% 
%         end
        feat(1,istants_of_common_presence)=cosangle(istants_of_common_presence);
        feat(2,istants_of_common_presence)=squeeze(sqrt(sum((skel1([1,3],3,istants_of_common_presence)-skel2([1,3],3,istants_of_common_presence)).^2)));
        feat(3,istants_of_common_presence)=squeeze(sqrt(sum((skel1([1,3],4,istants_of_common_presence)-skel2([1,3],3,istants_of_common_presence)).^2)));
        feat(4,istants_of_common_presence)=squeeze(sqrt(sum((skel1([1,3],6,istants_of_common_presence)-skel2([1,3],3,istants_of_common_presence)).^2)));
        feat(5,istants_of_common_presence)=squeeze(sqrt(sum((skel1([1,3],3,istants_of_common_presence)-skel2([1,3],4,istants_of_common_presence)).^2)));
        feat(6,istants_of_common_presence)=squeeze(sqrt(sum((skel1([1,3],3,istants_of_common_presence)-skel2([1,3],6,istants_of_common_presence)).^2)));
        Uorientdist{i,1}=feat;
        Fovdist{i,1}=dot(squeeze(skel2([1,2],3,:)),n1)-d1;
        Fovdist{i,2}=dot(squeeze(skel1([1,2],3,:)),n2)-d2;
        Fpdist{i,1}=dot(X,n1)-d1;
        Fpdist{i,2}=dot(X,n2)-d2;
        
    end


end



function [Denerg,MovingPixels,Variation] = Depthenerg(Depthclips,fps,w,th,dispth)
    if(nargin<4)
       th=4;
       dispth=0.1;
    end
    Denerg=cell(length(Depthclips.action),1);
    Variation=cell(length(Depthclips.action),1);
    MovingPixels=cell(length(Depthclips.action),1);
    for c=1:length(Depthclips.action)
        singleclip=cat(3,Depthclips.action{c}{:});
        Thresh=zeros(size(singleclip));
        Thresh(singleclip<th)=singleclip(singleclip<th);
        Delayed=circshift(Thresh,w,3);
        Delayed(:,:,1:w)=zeros(size(Delayed(:,:,1:w)));
        Variation{c}=abs(Delayed-Thresh);
        Denerg{c}=squeeze(sum(sum(abs(Variation{c}),2),1));
        MovingPixels{c}=squeeze(sum(sum(Variation{c}>dispth&Variation{c}<th)));
        
    end
end




function [ Venerg ] = energy_detailed(Sessions,fps,w)
    %fps = 15;
    %w = 5;
    

    numjoints=25;
    numdimensions=7;
    numsubjects=2;
    pos_columns=(reshape([1:7:175;(1:7:174)+1;(1:7:173)+2],1,length([1:7:175,(1:7:174)+1,(1:7:173)+2])));

    %noise_model=@(X)normpdf(X,0,0.01)*0.001;
    
    [numclips]=size(Sessions,1);
    Venerg=cell(numclips,numsubjects);
    for i=1:numclips

       lenclip=length(Sessions{i,2});
       skel1=nan(numdimensions,numjoints,lenclip);
       skel2=nan(numdimensions,numjoints,lenclip);
       skel1(:,:,Sessions{i,3}(1).temporal_presence)=reshape(Sessions{i,3}(1).data',numdimensions,numjoints,length(Sessions{i,3}(1).data));
       skel2(:,:,Sessions{i,3}(2).temporal_presence)=reshape(Sessions{i,3}(2).data',numdimensions,numjoints,length(Sessions{i,3}(2).data));
       skel1=skel1(1:3,:,:);
       skel2=skel2(1:3,:,:);
%        feat1=zeros(size(skel1,2),size(skel1,3));
%        feat2=zeros(size(skel2,2),size(skel2,3));
       feat1=zeros(size(skel1,2),lenclip);
       feat2=zeros(size(skel2,2),lenclip);
       istants_of_common_presence=find(Sessions{i,3}(1).temporal_presence&Sessions{i,3}(2).temporal_presence);
        
        for t=min(istants_of_common_presence):w+min(istants_of_common_presence)
            vel1=skel1(:,:,t)-repmat(skel1(:,3,t),1,numjoints);
            vel2=skel2(:,:,t)-repmat(skel2(:,3,t),1,numjoints);
            vel1(:,3)=skel1(:,3,t);
            vel2(:,3)=skel2(:,3,t);
            %vel1=vel1*w/fps;
            %vel2=vel2*w/fps;
            feat1(:,t)=sqrt(sum(vel1.^2));
            feat2(:,t)=sqrt(sum(vel2.^2));
            %feat1(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel1',ones(15,1),3))';
            %feat2(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel2',ones(15,1),3))';

        end
        for t=w+min(istants_of_common_presence):max(istants_of_common_presence)
            %vel1=skel1(:,:,t)-skel1(:,:,t-w);
            vel1=skel1(:,:,t)-repmat(skel1(:,3,t),1,numjoints)-skel1(:,:,t-w)+repmat(skel1(:,3,t-w),1,numjoints);
            vel1(:,3)=skel1(:,3,t)-skel1(:,3,t-w);
            %vel2=skel2(:,:,t)-skel2(:,:,t-w);
            vel2=skel2(:,:,t)-repmat(skel2(:,3,t),1,numjoints)-skel2(:,:,t-w)+repmat(skel2(:,3,t-w),1,numjoints);
            vel2(:,3)=skel2(:,3,t)-skel2(:,3,t-w);
            %vel1=vel1*w/fps;
            %vel2=vel2*w/fps;
            feat1(:,t)=sqrt(sum(vel1.^2));
            feat2(:,t)=sqrt(sum(vel2.^2));
            
%             feat1(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel1',ones(15,1),3))';
%             feat2(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel2',ones(15,1),3))';

        end
        Venerg{i,1}=feat1;
        Venerg{i,2}=feat2;
    end
end




function [ Ventr ] = vel_entropy_detailed(Sessions,fps,w)
%fps = 15;
%w = 5;
% pos_columns=(reshape([1:6:90;(1:6:89)+1;(1:6:88)+2],1,length([1:6:90,(1:6:89)+1,(1:6:88)+2])));

% noise_model=@(X)normpdf(X,0,0.01)*0.001;

numjoints=25;
numdimensions=7;
numsubjects=2;

[numclips]=size(Sessions,1);
Ventr=cell(numclips,numsubjects);
for i=1:numclips
    
   lenclip=length(Sessions{i,2});
   skel1=nan(numdimensions,numjoints,lenclip);
   skel2=nan(numdimensions,numjoints,lenclip);
   skel1(:,:,Sessions{i,3}(1).temporal_presence)=reshape(Sessions{i,3}(1).data',numdimensions,numjoints,length(Sessions{i,3}(1).data));
   skel2(:,:,Sessions{i,3}(2).temporal_presence)=reshape(Sessions{i,3}(2).data',numdimensions,numjoints,length(Sessions{i,3}(2).data));
   skel1=skel1(1:3,:,:);
   skel2=skel2(1:3,:,:);
   
   feat1=zeros(size(skel1,2),size(skel1,3));
   feat2=zeros(size(skel2,2),size(skel2,3));
   istants_of_common_presence=find(Sessions{i,3}(1).temporal_presence&Sessions{i,3}(2).temporal_presence);
        
   for t=min(istants_of_common_presence):w+min(istants_of_common_presence)
        vel1=skel1(:,:,t);
        vel2=skel2(:,:,t);
        feat1(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel1',ones(numjoints,1),3))';
        feat2(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel2',ones(numjoints,1),3))';
        
    end
    for t=w+min(istants_of_common_presence):max(istants_of_common_presence)
        %vel1=skel1(:,:,t)-skel1(:,:,t-w);
        vel1=skel1(:,:,t)-repmat(skel1(:,3,t),1,numjoints)-skel1(:,:,t-w)+repmat(skel1(:,3,t-w),1,numjoints);
        %vel2=skel2(:,:,t)-skel2(:,:,t-w);
        vel2=skel2(:,:,t)-repmat(skel2(:,3,t),1,numjoints)-skel2(:,:,t-w)+repmat(skel2(:,3,t-w),1,numjoints);
        feat1(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel1',ones(numjoints,1),3))';
        feat2(:,t)=cellfun(@(X)wentropy(X,'log energy'),mat2cell(vel2',ones(numjoints,1),3))';
        
    end
    Ventr{i,1}=feat1;
    Ventr{i,2}=feat2;
end
%


end





function [ Ventr ] = vel_entropy( Sessions, fps, w)
%VEL_ENTROPY Converts data points of the skeleton in velocity entropy


velxyz = struct('activity', {});
%fps = 15;
%w = 5;

pos_columns=(reshape([1:7:175;(1:7:174)+1;(1:7:173)+2],1,length([1:7:175,(1:7:174)+1,(1:7:173)+2])));

[numclips]=size(Sessions,1);
numsubjects=2;
Ventr=cell(numclips,3);

    for i=1:numclips
       clip=Sessions{i,3};
       lenclip=length(Sessions{i,2});
       
       skel1=nan(lenclip,length(pos_columns));
       skel1(Sessions{i,3}(1).temporal_presence,:)=Sessions{i,3}(1).data(Sessions{i,3}(1).temporal_presence,pos_columns);
       skel2=nan(lenclip,length(pos_columns));
       skel2(Sessions{i,3}(2).temporal_presence,:)=Sessions{i,3}(2).data(Sessions{i,3}(2).temporal_presence,pos_columns);
       
       vele=nan(lenclip,1);
       vele1=nan(lenclip,1);
       vele2=nan(lenclip,1);
       istants_of_common_presence=find(Sessions{i,3}(1).temporal_presence&Sessions{i,3}(2).temporal_presence);
        for t=(w+min(istants_of_common_presence)):max(istants_of_common_presence)
%             disp(t)
            vele(t)=wentropy([skel1(t,:),skel2(t,:)]-[skel1(t-w,:),skel2(t-w,:)],'log energy');
            vele1(t)=wentropy(skel1(t,:)-skel1(t-w,:),'log energy');
            vele2(t)=wentropy(skel2(t,:)-skel2(t-w,:),'log energy');
        end
        
        Ventr{i,1}=vele;
        Ventr{i,2}=vele1;
        Ventr{i,3}=vele2;
    end
    
    
    
end

