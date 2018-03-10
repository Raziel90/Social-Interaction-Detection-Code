clear
clc


k2_init
class_seg_path='/home/ccoppola/Social_Detection_Recognition';
file_saving_path='/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/Useful_files';




[Class_Annotation,classnames,filename]=k2_importclassAnnotationfile([class_seg_path,'/session_class_1/Annotation.csv'],1, 11);
skeleton_class_bags=k2_skeleton_text_to_matlab( [class_seg_path,'/Output'] );
% for k=1:size(skeleton_class_bags)
% 
% skel1=nan(length(skeleton_class_bags{k,2}),7*25);
% skel2=nan(length(skeleton_class_bags{k,2}),7*25);
% mutual_presence=skeleton_class_bags{k,3}(1).temporal_presence&skeleton_class_bags{k,3}(2).temporal_presence;
% skel1(skeleton_class_bags{k,3}(1).temporal_presence,:)=skeleton_class_bags{3}(1).data;
% skel2(skeleton_class_bags{k,3}(2).temporal_presence,:)=skeleton_class_bags{3}(2).data;
% skeleton_class_bags{k,3}(1).data=skel1(skeleton_class_bags{3}(1).temporal_presence&skeleton_class_bags{3}(2).temporal_presence,:);
% skeleton_class_bags{k,3}(2).data=skel2(skeleton_class_bags{3}(1).temporal_presence&skeleton_class_bags{3}(2).temporal_presence,:);
% skeleton_class_bags{k,2}=skeleton_class_bags{k,2}(skeleton_class_bags{k,3}(1).temporal_presence&skeleton_class_bags{k,3}(2).temporal_presence);
% skeleton_bags{k,3}(1).temporal_presence()
% end
Class_Truth=zeros(size(skeleton_class_bags{2}));
for h=1:length(classnames)
    
    for k=1:2:(sum(~isnan(Class_Annotation(:,h))))
        Class_Truth(skeleton_class_bags{2}<Class_Annotation(k+1,h)&skeleton_class_bags{2}>=Class_Annotation(k,h))=h;
        
    end
    
    
end


clear h k skel1 skel2








%Temporary variable
wrap_Gtruth.segmentation_user1=double(Class_Truth>0)';
wrap_Gtruth.segmentation_user2=double(Class_Truth>0)';
Class_Truth=Class_Truth(skeleton_class_bags{1,3}(1).temporal_presence&skeleton_class_bags{1,3}(2).temporal_presence);


[input_seg_test_feat,input_seg_test_lab]=k2_compute_features(skeleton_class_bags,wrap_Gtruth);
clear h k skel1 skel2 wrap_Gtruth

[input_seg_tr_feat,input_seg_tr_lab]=k2_compute_features(skeleton_bags,G_Truth);


%%MIN-MAX NORMAILZATION
max_tr=max(cell2mat(cellfun(@(X)max(X,[],2),{input_seg_tr_feat{:}},'UniformOutput',false)),[],2)';
min_tr=min(cell2mat(cellfun(@(X)min(X,[],2),{input_seg_tr_feat{:}},'UniformOutput',false)),[],2)';

for k=1:size(input_seg_tr_feat,1)
    input_seg_tr_feat{k}(max_tr~=min_tr,:)=(input_seg_tr_feat{k}(max_tr~=min_tr,:)-min_tr(max_tr~=min_tr)')./(max_tr(max_tr~=min_tr)-min_tr(max_tr~=min_tr))';
end
for k=1:size(input_seg_test_feat,1)
    input_seg_test_feat{k}(max_tr~=min_tr,:)=(input_seg_test_feat{k}(max_tr~=min_tr,:)-min_tr(max_tr~=min_tr)')./(max_tr(max_tr~=min_tr)-min_tr(max_tr~=min_tr))';
end


%%DOWNSAMPLING
tr_X=horzcat(input_seg_tr_feat{:})';
tr_lab=horzcat(input_seg_tr_lab{:})';
te_X=horzcat(input_seg_test_feat{:})';
te_lab=horzcat(input_seg_test_lab{:})';

% if 1==1
%    downsampling=1:3:size(tr_X,1);
%     tr_X=tr_X(downsampling,:);
%     tr_lab=tr_lab(downsampling);
% end


%%SEGMENTATION
gamma=1;
cost=1;
kernel_code=0;
[hmms.TRANS,hmms.INIT]=transmat_train_observed({input_seg_tr_lab{:}}',2);
[Btr, Bte, cm_tr_svm, cm_te_svm, hmms.EMISS] = SVM_tr_te(tr_X,te_X,tr_lab,te_lab,kernel_code,cost,gamma);
% Accuracy = 87.4867% (41103/46982) (classification)
% Accuracy = 79.5812% (912/1146) (classification)
seg_lab=viterbi_path(hmms.INIT, hmms.TRANS,Bte');
seg_results=seg_lab-1;
[ segm_stats{1},segm_stats{2},h] = classification_stats( input_seg_test_lab{:}-1,horzcat(seg_results),1,{'individual','social'});

ideal_start_stop=find(abs(diff(input_seg_test_lab{:}-1))~=0);
if mod(size(ideal_start_stop,1),2)~=0&&te_lab(end)==1
    ideal_start_stop=[ideal_start_stop,length(te_lab)];
elseif mod(size(ideal_start_stop,1),2)~=0&&te_lab(1)==1
    ideal_start_stop=[1,ideal_start_stop];
end
ideal_start_stop=reshape(ideal_start_stop,2,length(ideal_start_stop)/2);
ideal_start_stop(1,:)=ideal_start_stop(1,:)+1;
real_start_stop=find(abs(diff(seg_results))~=0);

if mod(size(real_start_stop,1),2)~=0&&seg_results(end)==1
    real_start_stop=[real_start_stop,length(seg_results)];
elseif mod(size(real_start_stop,1),2)~=0&&seg_results(1)==1
    real_start_stop=[1,real_start_stop];
end
real_start_stop=reshape(real_start_stop,2,length(real_start_stop)/2);
real_start_stop(1,:)=real_start_stop(1,:)+1;

idseg_social_priors_feat=repmat(struct('DistTT',[],'MinDistT1',[],'MinDistT2',[],'MaxDistT1',[],'MaxDistT2',[],'MinDist',[],'MaxDist',[]),size(ideal_start_stop,2),2);
realseg_social_priors_feat=repmat(struct('DistTT',[],'MinDistT1',[],'MinDistT2',[],'MaxDistT1',[],'MaxDistT2',[],'MinDist',[],'MaxDist',[]),size(real_start_stop,2),1);
idseg_test_social_feat=repmat(struct('activity',[]),size(ideal_start_stop,2),2);
realseg_test_social_feat=repmat(struct('activity',[]),size(real_start_stop,2),2);
IDTest_Sessions=repmat(struct('activity',[]),size(ideal_start_stop,2),2);
REALTest_Sessions=repmat(struct('activity',[]),size(real_start_stop,2),2);
%%CLASSIFICATION
load('Dataset_Social_features_priors_v2.mat')
%Downsample to 15fps
for p=1:size(Sessions,1)
    for a=1:size(Sessions(p,1).action,1)
        Sessions(p,1).action{a}=Sessions(p,1).action{a}(1:2:end,:);
        Sessions(p,2).action{a}=Sessions(p,2).action{a}(1:2:end,:);
    end
end

tr_social_feat=sk_social_features(10,8,Sessions(:,1),Sessions(:,2));
tr_social_priors_feat=sk_social_priors_features(10,8,Sessions(:,1),Sessions(:,2));
[priors_dist,multivariate_dist]=sk_priors_distr_learn(8,tr_social_priors_feat,4);

[Test_Sessions]=k2_data_retrocompatibility(skeleton_class_bags);
test_social_feat=sk_social_features(1,1,Test_Sessions(1),Test_Sessions(2));
test_social_priors_feat=sk_social_priors_features(1,1,Test_Sessions(1),Test_Sessions(2));
[multivariate_unseg_priors,multivariate_unseg_posteriors] = sk_multivariate_social_priors(1, 8, multivariate_dist,test_social_priors_feat,1);
[unseg_priors,unseg_posteriors] = sk_social_priors(1, 1, priors_dist,test_social_priors_feat);
% [priors_dist,multivariate_dist]=sk_priors_distr_learn(1,test_social_priors_feat,4);

%%APPLY SEGMENTATION TO FEATURES
for k=1:size(ideal_start_stop,2)
    idseg_social_priors_feat(k).DistTT={test_social_priors_feat.DistTT{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MinDistT1={test_social_priors_feat.MinDistT1{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MinDistT2={test_social_priors_feat.MinDistT2{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MaxDistT1={test_social_priors_feat.MaxDistT1{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MaxDistT2={test_social_priors_feat.MaxDistT2{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MinDist={test_social_priors_feat.MinDist{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_social_priors_feat(k).MaxDist={test_social_priors_feat.MaxDist{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    idseg_test_social_feat(k).activity={test_social_feat(1).activity{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    IDTest_Sessions(k,1).action={Test_Sessions(1).action{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
    IDTest_Sessions(k,2).action={Test_Sessions(2).action{1}(ideal_start_stop(1,k):ideal_start_stop(2,k),:)};
end
[idseg_priors,idseg_posteriors]=sk_social_priors(size(ideal_start_stop,2),1,priors_dist,idseg_social_priors_feat);
[multivariate_idseg_priors,multivariate_idseg_posteriors] = sk_multivariate_social_priors(size(ideal_start_stop,2), 8, multivariate_dist,idseg_social_priors_feat,1);

for k=1:size(real_start_stop,2)
    realseg_social_priors_feat(k).DistTT={test_social_priors_feat.DistTT{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MinDistT1={test_social_priors_feat.MinDistT1{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MinDistT2={test_social_priors_feat.MinDistT2{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MaxDistT1={test_social_priors_feat.MaxDistT1{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MaxDistT2={test_social_priors_feat.MaxDistT2{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MinDist={test_social_priors_feat.MinDist{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_social_priors_feat(k).MaxDist={test_social_priors_feat.MaxDist{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    realseg_test_social_feat(k).activity={test_social_feat(1).activity{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    REALTest_Sessions(k,1).action={Test_Sessions(1).action{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
    REALTest_Sessions(k,2).action={Test_Sessions(2).action{1}(real_start_stop(1,k):real_start_stop(2,k),:)};
end

[realseg_priors,realseg_posteriors]=sk_social_priors(size(real_start_stop,2),1,priors_dist,realseg_social_priors_feat);
[multivariate_realseg_priors,multivariate_realseg_posteriors] = sk_multivariate_social_priors(size(real_start_stop,2), 8, multivariate_dist,realseg_social_priors_feat,1);


% Sessions=Test_Sessions;
% social_features=test_social_feat;
% social_priors_feat=test_social_priors_feat;
% multivariate_priors=multivariate_unseg_priors;
% multivariate_posteriors=multivariate_unseg_posteriors;
% priors=unseg_priors;
% posteriors=unseg_posteriors;
% save([file_saving_path,'/k2_unsegmented_new.mat'],'Sessions','Class_Truth','tr_social_feat','tr_social_priors_feat','priors_dist','multivariate_dist','social_features','test_social_priors_feat','priors','posteriors','multivariate_priors','multivariate_posteriors');
save([file_saving_path,'/k2_unsegmented.mat'],'Sessions','Test_Sessions','Class_Truth','tr_social_feat','tr_social_priors_feat','priors_dist','multivariate_dist','test_social_feat','test_social_priors_feat','unseg_priors','unseg_posteriors','multivariate_unseg_priors','multivariate_unseg_posteriors');
save([file_saving_path,'/k2_idsegmented.mat'],'Sessions','IDTest_Sessions','Class_Truth','tr_social_feat','tr_social_priors_feat','ideal_start_stop','priors_dist','multivariate_dist','test_social_feat','test_social_priors_feat','idseg_priors','idseg_posteriors','multivariate_realseg_priors','multivariate_idseg_posteriors');
save([file_saving_path,'/k2_realsegmented.mat'],'Sessions','REALTest_Sessions','Class_Truth','tr_social_feat','tr_social_priors_feat','real_start_stop','priors_dist','multivariate_dist','test_social_feat','test_social_priors_feat','realseg_priors','realseg_posteriors','multivariate_idseg_priors','multivariate_realseg_posteriors');









