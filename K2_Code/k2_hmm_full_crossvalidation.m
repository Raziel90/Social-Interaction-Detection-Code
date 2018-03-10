%% hard parameters
numdimensions=7;
numjoints=25;


%% soft parameters
if ~exist('select_parameters')
distance_sectors=[0.15,0.46,0.76,1.22,2.1,3.7,7.6];
% distance_sectors=[0.45,1.2,3.7,7.6]; %chosen according to proxemics theory
mode=5;
% distance_sectors=[];
QTC_threshold=0.3/15;
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
% feat_sel=logical([1,0,0,0,0,0]);
feat_str=num2str(feat_sel,[repmat('%i_',1,length(feat_sel)-1),'%i']);

covariance_applied=feat_sel(6);
eig_mode=0;

end


%% init
gt=cell(10,1);
hmms=cell(10,1);
results=cell(10,1);
input_sequences=cell(10,1);
nan_values=cell(10,1);
test_probabilities=cell(10,1);
training_probabilities=cell(10,1);
super_dist=cell(10,1);
real_dist=cell(10,1);
qtc_windowed_strings=cell(10,1);
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
    nan_values1=squeeze(isnan(skel1(1,1,:)));
    nan_values2=squeeze(isnan(skel2(1,1,:)));
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

%% Leaveoneout cv
for j=1:size(skeleton_bags,1)
    
    testing_clip=all_clips(j);
    training_clip=all_clips(all_clips~=j);
    

    if mode >=3
        % Normalization
        max_tr=max(cell2mat(cellfun(@(X)max(X,[],2),{input_sequences{training_clip}},'UniformOutput',false)),[],2)';
        min_tr=min(cell2mat(cellfun(@(X)min(X,[],2),{input_sequences{training_clip}},'UniformOutput',false)),[],2)';
        for k=1:length(all_clips)
            input_sequences{k}(max_tr~=min_tr,:)=(input_sequences{k}(max_tr~=min_tr,:)-min_tr(max_tr~=min_tr)')./(max_tr(max_tr~=min_tr)-min_tr(max_tr~=min_tr))';
        end
    end
    
    tr_X=horzcat(input_sequences{training_clip})';
    tr_lab=horzcat(gt{training_clip})';
    te_X=input_sequences{testing_clip}';
    te_lab=gt{testing_clip}';
%     for k=1:size(skeleton_bags,1)

        if pcait==1
            [pc,score,latent,~,explained] = pca(tr_X);
            [~,lastcomponent]=min(abs(cumsum(explained)-95));
            %[tr_X,mu,sigma]=zscore(tr_X);
            tr_X=tr_X*pc(:,1:lastcomponent);
%             te_X=(te_X-mu./sigma)*pc(:,1:lastcomponent);
            te_X=te_X*pc(:,1:lastcomponent);
        end
        if subsample==1
           downsampling=1:3:size(tr_X,1);
            tr_X=tr_X(downsampling,:);
            tr_lab=tr_lab(downsampling);
        end


        if(mode==1)
            [hmms{j}.TRANS,hmms{j}.EMISS]=hmmestimate(tr_X,tr_lab,'Pseudoemissions',ones(2,num_of_possible_inputs)*1e-10,'Pseudotransitions',ones(2,2)*1e-10);
            hmms{j}.TRANS=[0,[0.5,0.5];zeros(size(hmms{j}.TRANS,1),1),hmms{j}.TRANS];
            hmms{j}.EMISS=[zeros(1,size(hmms{j}.EMISS,2));hmms{j}.EMISS];
            classstr='MultinomialHMM';
        elseif mode==0
           [hmms{j}.TRANS,hmms{j}.EMISS]=hmmtrain({input_sequences{training_clip}},[0,[0.5,0.5];zeros(2,1),rand(2)],[zeros(1,num_of_possible_inputs);rand(2,num_of_possible_inputs)],'ALGORITHM','BaumWelch','Pseudoemissions',ones(3,num_of_possible_inputs)*1e-10,'Pseudotransitions',ones(3,3)*1e-10);
           classstr='MultinomialHMM_EM';
        elseif mode==2
            [hmms{j}.TRANS,hmms{j}.EMISS]=hmmestimate(tr_X,tr_lab,'Pseudoemissions',ones(2,num_of_possible_inputs)*1e-10,'Pseudotransitions',ones(2,2)*1e-10);
            [hmms{j}.TRANS,hmms{j}.EMISS]=hmmtrain({input_sequences{training_clip}},hmms{j}.TRANS,hmms{j}.EMISS,'ALGORITHM','BaumWelch','Pseudoemissions',ones(2,num_of_possible_inputs)*1e-10,'Pseudotransitions',ones(2,2)*1e-10);
            hmms{j}.TRANS=[0,[0.5,0.5];zeros(size(hmms{j}.TRANS,1),1),hmms{j}.TRANS];
            hmms{j}.EMISS=[zeros(1,size(hmms{j}.EMISS,2));hmms{j}.EMISS];
            classstr='MultinomialHMM_EMenforced';
        elseif mode==3
            classstr='MGaussianHMM';
            [hmms{j}.INIT, hmms{j}.TRANS, hmms{j}.MU, hmms{j}.SIGMA]=gausshmm_train_observed({input_sequences{training_clip}}',{gt{training_clip}}',2,'cov_type','full');
            hmms{j}.MIX=ones(2,1);
        elseif mode==4
            classstr='GMM_HMM';
            [hmms{j}.INIT, hmms{j}.TRANS, hmms{j}.MU, hmms{j}.SIGMA]=mhmmParzen_train_observed({input_sequences{training_clip}}',{gt{training_clip}}',2,3);
            mixgauss = mixgauss_classifier_train(horzcat(input_sequences{training_clip}), horzcat(gt{training_clip})-1, 5,'cov_type','full','max_iter', 100);
            hmms{j}.MIX=[mixgauss.neg.prior';mixgauss.pos.prior'];
            hmms{j}.MU=permute(cat(3,mixgauss.neg.mu,mixgauss.pos.mu),[1,3,2]);
            hmms{j}.SIGMA=permute(cat(4,mixgauss.neg.Sigma,mixgauss.pos.Sigma),[1,2,4,3]);
        elseif mode==5
            classstr='SVM_HMM';
            
            [hmms{j}.TRANS,hmms{j}.INIT]=transmat_train_observed({gt{training_clip}}',2);
            [Btr, Bte, cm_tr_svm, cm_te_svm, hmms{j}.EMISS] = SVM_tr_te(tr_X,te_X,tr_lab,te_lab,kernel_code,cost,gamma);
        elseif mode==6
            classstr='RandomF_HMM';
            [hmms{j}.TRANS,hmms{j}.INIT]=transmat_train_observed({gt{training_clip}}',2);
            opts.depth= 11;
            opts.numTrees= 300; %(but more is _ALWAYS_ better, monotonically, no exceptions)
            opts.numSplits= 10;
            opts.classifierID= 1;
            opts.classifierCommitFirst= true;
            model=forestTrain(tr_X,tr_lab,opts);
            [~,Bte]=forestTest(model,te_X);
            [~,Btr]=forestTest(model,tr_X);
            %[Btr, Bte, cm_tr_svm, cm_te_svm, hmms{j}.EMISS] = SVMclass(training,input_sequences{testing_clip}',labels,gt{testing_clip}');
        end
%     end
    if mode<3&&mode>0
        [pstates,logpseq]=hmmdecode(input_sequences{testing_clip},hmms{j}.TRANS,hmms{j}.EMISS);
        [~,s]=max(pstates,[],1);
        results{j}=s-2;
    elseif mode<5
        [B, B2] = mixgauss_prob(te_X', hmms{j}.MU, hmms{j}.SIGMA, hmms{j}.MIX);   %, unit_norm)
        s=viterbi_path(hmms{j}.INIT, hmms{j}.TRANS,B);
        results{j}=s-1;
%         [loglik,errors]=mhmm_logprob(input_sequences{testing_clip},hmms{j}.INIT,hmms{j}.TRANS,hmms{j}.MU,hmms{j}.SIGMA,hmms{j}.MIX);
    elseif mode>=5
        test_probabilities{j}=Bte;
        training_probabilities{j}=Btr;
        s=viterbi_path(hmms{j}.INIT, hmms{j}.TRANS,Bte');
        results{j}=s-1;
    end
disp({['cross-validation-iter: ',num2str(j)];['Classifier: ',classstr, ' datasize: ',num2str(size(tr_X))]})
end

closed_results=cell(size(results));

for k=1:length(results)
    closed_results{k}=imerode(imdilate(results{k},strel('line',closing_amount,1)),strel('line',closing_amount,1));
end

%% plot results
correct=0;
total=0;
for k=1:size(skeleton_bags,1)
    correct=correct+sum(results{k}==(gt{k}-1));
total=total+length(results{k});
end

per_clip_accuracy=cell(size(results,1),2);
for k=1:size(skeleton_bags,1)
%     per_clip_accuracy{k}=sum(results{k}==(gt{k}-1))/length(results{k});
    [ per_clip_accuracy{k,1},per_clip_accuracy{k,2}] = classification_stats( gt{k}-1,results{k},0);
end

folname=[classstr,'_K',num2str(kernel_code),'_G',num2str(round(log2(gamma))),'_C',num2str(round(log2(cost))),'/'];
succ=mkdir([outputfolder,'/k2_figures/k2_',classstr,'_feats',feat_str,pcastr(pcait&true(1,length(pcastr))),'/'],folname);
for k=1:size(skeleton_bags,1)
    
    h=figure(k); plot(skeleton_bags{k,2}(~nan_values{k}),results{k},'b'),hold on,plot(skeleton_bags{k,2}(~nan_values{k}),(gt{k}-1)'*1.1,'r');
    title(['clip-',num2str(k)])
    if saveit==1
        
        
        savefig(h,[outputfolder,'/k2_figures/k2_',classstr,'_feats',feat_str,pcastr(pcait&true(1,length(pcastr))),'/',folname,['qtcsize',num2str(qtc_size),'_prox','_angle_',classstr],'result_clip',num2str(k)]);
        print([outputfolder,'/k2_figures/k2_',classstr,'_feats',feat_str,pcastr(pcait&true(1,length(pcastr))),'/',folname,['qtcsize',num2str(qtc_size),'_prox','_angle_',classstr],'result_clip',num2str(k)], '-dpng', '-r300');
    end
end
close all;
% for k=1:size(skeleton_bags,1)
%     
%     h=figure(k); plot(skeleton_bags{k,2}(~nan_values{k}),closed_results{k},'b'),hold on,plot(skeleton_bags{k,2}(~nan_values{k}),(gt{k}-1)'*1.1,'r');
%     title(['clip-',num2str(k)])
%     if saveit==1
%     %     saveitfig(h,['/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/figures/k2_hmm(qtc+prox)/closed_result_clip',num2str(k)]);
% %         print(['/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/figures/k2_hmm(qtc+prox)/closed_result_clip',num2str(k)], '-dpng', '-r300');
%     end
% end
close all;
[ class_stats{1},class_stats{2},h] = classification_stats( horzcat(gt{:})-1,horzcat(results{:}),1,{'individual','social'});
% h=figure; plotconfusion(horzcat(results{:}),horzcat(gt{:})-1);
if saveit==1
    savefig(h,[outputfolder,'/k2_figures/k2_',classstr,'_feats',feat_str,pcastr(pcait&true(1,length(pcastr))),'/',folname,['qtcsize',num2str(qtc_size),'_prox','_angle_mode',classstr],'confusion_mat']);
    print([outputfolder,'/k2_figures/k2_',classstr,'_feats',feat_str,pcastr(pcait&true(1,length(pcastr))),'/',folname,['qtcsize',num2str(qtc_size),'_prox','_angle_mode',classstr],'confusion_mat'], '-dpng', '-r300');
end
close all;

clear s j k skel1 skel2 labels training nan_values1 nan_values2 total training_clip testing_clip lenclip input_sequences correct cm_te_svm cm_tr_svm B B2 Btr Bte folname succ