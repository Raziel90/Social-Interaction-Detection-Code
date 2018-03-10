function [input_sequences,gt]=k2_compute_features(skeleton_bags,G_Truth,feat_sel)

%% hard parameters
numdimensions=7;
numjoints=25;


%% soft parameters
if ~exist('select_parameters')
distance_sectors=[0.15,0.46,0.76,1.22,2.1,3.7,7.6];
% distance_sectors=[0.45,1.2,3.7,7.6]; %chosen according to proxemics theory
mode=5;
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

gamma=1;
cost=1;
kernel_code=0;

%downsampling the trainig
subsample=1;
pcait=0;pcastr='withpca';
%feature selection
feat_sel=logical([1,1,1,1,1,1]);
feat_str=num2str(feat_sel,[repmat('%i_',1,length(feat_sel)-1),'%i']);

covariance_applied=feat_sel(6);
eig_mode=0;

end


%% init
gt=cell(size(skeleton_bags,1),1);
hmms=cell(size(skeleton_bags,1),1);
results=cell(size(skeleton_bags,1),1);
input_sequences=cell(size(skeleton_bags,1),1);
nan_values=cell(size(skeleton_bags,1),1);
test_probabilities=cell(size(skeleton_bags,1),1);
training_probabilities=cell(size(skeleton_bags,1),1);
super_dist=cell(size(skeleton_bags,1),1);
real_dist=cell(size(skeleton_bags,1),1);
classstr='';

%% Compute features
[qtc_sequences] = k2_convert_qtc_all(skeleton_bags,QTC_threshold,2,mode<3);
[dist_sequences]=k2_discretize_distance(skeleton_bags, distance_sectors,hysteresis, 5);

num_of_possible_qtc_codes=prod(repmat(3,1,qtc_size));

num_of_possible_inputs=prod([length(distance_sectors)+1,num_of_possible_qtc_codes]);


qtc_hists=k2_QTC_moving_hist(skeleton_bags,15,eps,2,0);


%[orient12,orient21] = k2_compute_body_orientation(skeleton_bags,6,10,0);
[dir2centre1,dir2centre2,radius,orient12,orient21] = k2_overlap_transactional_segment(skeleton_bags,6,10,0);
for k=1:size(skeleton_bags,1)
    lenclip=length(skeleton_bags{k,2});
    skel1=nan(numdimensions,numjoints,lenclip);
    skel2=nan(numdimensions,numjoints,lenclip);
    skel1(:,:,skeleton_bags{k,3}(1).temporal_presence)=reshape(skeleton_bags{k,3}(1).data',numdimensions,numjoints,length(skeleton_bags{k,3}(1).data));
    skel2(:,:,skeleton_bags{k,3}(2).temporal_presence)=reshape(skeleton_bags{k,3}(2).data',numdimensions,numjoints,length(skeleton_bags{k,3}(2).data));
    Jdist=cellfun(@(X,Y)dist2(X',Y'),squeeze(num2cell(skel1([1,2],[3,5,6,10],:),[1,2])),squeeze(num2cell(skel2([1,2],[3,5,6,10],:),[1,2])),'UniformOutput',false);
    super_dist{k}=reshape(cat(3,Jdist{:}),4*4,length(Jdist));
    skel1=skel1([1,2],3,:);
    skel2=skel2([1,2],3,:);
%     nan_values1=squeeze(isnan(skel1(1,1,:)));
%     nan_values2=squeeze(isnan(skel2(1,1,:)));
    nan_values1=~skeleton_bags{k,3}(1).temporal_presence;
    nan_values2=~skeleton_bags{k,3}(2).temporal_presence;
    nan_values{k}=nan_values1|nan_values2;
    super_dist{k}=super_dist{k}(:,~nan_values{k});
    real_dist{k}=medfilt1(sqrt(sum(squeeze(skel1(:,1,~nan_values{k})-skel2(:,1,~nan_values{k})).^2)),med_filt_lev);
    super_dist{k}=medfilt1(super_dist{k},med_filt_lev,[],2);
    gt{k}=1+double(G_Truth(k).segmentation_user1(~nan_values{k})&G_Truth(k).segmentation_user2(~nan_values{k}))';
    qtc_sequences{k}=qtc_sequences{k}(~nan_values{k},1:min(qtc_size,size(qtc_sequences{k},2)));
	
    wsize=1;
    if wsize>1
        windowed_qtc=zeros(size(qtc_sequences{k},1),qtc_size*wsize);
        windowed_qtc(:,qtc_size*(wsize-1)+1:qtc_size*wsize)=qtc_sequences{k};
        for qtccat=(wsize-2):-1:0
            windowed_qtc(1:size(qtc_sequences{k},1)-wsize+qtccat+1,qtc_size*qtccat+1:qtc_size*(qtccat+1))=qtc_sequences{k}(wsize-qtccat:size(qtc_sequences{k},1),:);
        end
        qtc_sequences{k}=windowed_qtc;
        qtc_windowed_strings{k}=cell2mat(cellfun(@(Y)strrep(strrep(strrep(Y','-1','-'),'1','+'),' ','')',mat2cell(num2str(windowed_qtc')',size(num2str(windowed_qtc')',1),ones(1,size(windowed_qtc,2))),'UniformOutput',false));
%         qtc_windowed_strings{k}=cell2mat(cellfun(@(Y)strrep(strrep(strrep(Y','-1','-'),'1','+'),' ','')',mat2cell(char(33+windowed_qtc')',size(char(33+windowed_qtc')',1),ones(1,size(windowed_qtc,2))),'UniformOutput',false));
    end
    dist_sequences{k}=dist_sequences{k}(~nan_values{k});
    orient12{k}=medfilt1(orient12{k}(~nan_values{k}),med_filt_lev);
    orient21{k}=medfilt1(orient21{k}(~nan_values{k}),med_filt_lev);
    qtc_hists{k}=qtc_hists{k}(~nan_values{k},:);
    dir2centre1{k}=dir2centre1{k}(~nan_values{k});
    dir2centre2{k}=dir2centre2{k}(~nan_values{k});
    radius{k}=radius{k}(~nan_values{k});
    vel{k}=[0,diff(real_dist{k})];
%     gt{k}=1+double(G_Truth(k).segmentation_user1&G_Truth(k).segmentation_user2)';
    
end
%% select input
for k=1:size(skeleton_bags,1)
    if mode>0&&mode<3
    if qtc_size>0 && length(distance_sectors)>2
    	input_sequences{k}=qtc_sequences{k}'+3^4.*dist_sequences{k};
        filename_base=['qtcc_prox',num2str(length(distance_sectors)+1),'_'];
    elseif qtc_size==0 && length(distance_sectors)>2
        input_sequences{k}=dist_sequences{k}+1;
        filename_base=['prox',num2str(length(distance_sectors)+1),'_'];
    elseif qtc_size>0 && length(distance_sectors)<2
        input_sequences{k}=qtc_sequences{k}';
        filename_base=['qtcc_'];
    end
    elseif mode>=3

        
%         input_sequences{k}=[qtc_sequences{k};medfilt1(sqrt(sum(squeeze(skel1(:,1,~nan_values)-skel2(:,1,~nan_values)).^2)),20)];
        input_sequences{k}=[];
        
%         input_sequences{k}=[input_sequences{k};qtc_sequences{k}];
        input_sequences{k}=[input_sequences{k};qtc_sequences{k}(repmat(feat_sel(1),1,max(size(qtc_sequences{k}))),:)'];
%         input_sequences{k}=[input_sequences{k};real_dist{k}];
        %input_sequences{k}=[input_sequences{k};real_dist{k}(:,repmat(feat_sel(2),1,max(size(real_dist{k}))))];
        input_sequences{k}=[input_sequences{k};super_dist{k}(:,repmat(feat_sel(2),1,max(size(super_dist{k}))))];
%         input_sequences{k}=[input_sequences{k};orient12{k};orient21{k}];
        input_sequences{k}=[input_sequences{k};orient12{k}(:,repmat(feat_sel(3),1,max(size(orient12{k}))));orient21{k}(:,repmat(feat_sel(3),1,max(size(orient21{k}))))];
%         input_sequences{k}=[input_sequences{k};qtc_hists{k}];
        input_sequences{k}=[input_sequences{k};qtc_hists{k}(repmat(feat_sel(4),1,max(size(qtc_hists{k}))),:)'];
%         input_sequences{k}=[input_sequences{k};dir2centre1{k};dir2centre2{k};radius{k}];
        input_sequences{k}=[input_sequences{k};dir2centre1{k}(:,repmat(feat_sel(5),1,max(size(dir2centre1{k}))));dir2centre2{k}(:,repmat(feat_sel(5),1,max(size(dir2centre2{k}))));radius{k}(:,repmat(feat_sel(5),1,max(size(radius{k}))))];
        
        
    end
    
    
end

if covariance_applied==1&&feat_sel(6)
    %cov_seq=k2_data_moving_cov(input_sequences, 15, [2,3],eig_mode);
    cov_seq=k2_data_moving_cov(cellfun(@(X,Y)vertcat(X,Y),orient12,orient21,'UniformOutput',false), 15, [1,2],eig_mode);
    for k=1:length(cov_seq)    
        input_sequences{k}=[input_sequences{k};cov_seq{k}];
    end
end
end
