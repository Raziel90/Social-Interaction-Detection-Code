clear 
clc

        
k2_init



%STATIC CPARAMS

%num of tested features
feat_num=6;

cost_choices=2.^[-5:2:5];
kernel_code=0;
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
mode=5;
subsample=1;
%kernel_code,cost,gamma


select_parameters=1;



distance_sectors=[0.15,0.46,0.76,1.22,2.1,3.7,7.6];
% distance_sectors=[0.45,1.2,3.7,7.6]; %chosen according to proxemics theory

% distance_sectors=[];
% QTC_threshold=0.3/15;
QTC_threshold=eps;
all_clips=1:size(skeleton_bags,1);
qtc_size=4;
closing_amount=45;
hysteresis=[0.0,0.5];
saveit=1;
med_filt_lev=20;
filename_base='';
covariance_applied=1;
eig_mode=0;

pcait=0;pcastr='withpca';


comparison=cell(feat_num,1);

for q=1:length(cost_choices)%feat_num
    gamma=0;
%     cost=cost_choices(q);
    cost=1;
    feat_sel=logical(zeros(1,feat_num));
    feat_sel(q)=or(feat_sel(q),1);
    feat_str=num2str(feat_sel,[repmat('%i_',1,length(feat_sel)-1),'%i']);
    k2_hmm_full_crossvalidation
    comparison{q}=class_stats;
    save([basepath,'Social-Activities-Code/k2_results/',classstr,'configuration_K',num2str(kernel_code),'_c',num2str(round(log2(cost))),'_feat',feat_str,'.mat'],'hmms','classstr','results','per_clip_accuracy','kernel_code','cost','gamma','mode','gt','G_Truth','closed_results','closing_amount','nan_values','test_probabilities','training_probabilities','class_stats');
    clear('hmms','classstr','results','per_clip_accuracy','cost','gamma','gt','closed_results','nan_values','test_probabilities','training_probabilities','class_stats');
end
comparison=cat(1,comparison{:});
clear select_parameters feat_num