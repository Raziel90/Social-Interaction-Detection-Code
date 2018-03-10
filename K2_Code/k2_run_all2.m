


function k2_run_all()
% addpath(genpath('~/develop/Social-Behaviour/'))

    tic
    clc;
    clear;
    num_persons = 10;
	num_tests = 1;
	numact = 8;
    %sk_matfile = 'Dataset_Social_features_priors_v2.mat';
    sk_matfile = '2016_Social_Dataset_downsampled15fps.mat';
    sk_test_matfile='k2_unsegmented.mat';
    % indexes 1-'Bathroom'; 2-'Bedroom'; 3-'Kitchen'; 4-'Livingroom';
    % 5-'Office'
    %rooms_desc = ['Bathroom  '; 'Bedroom   '; 'Kitchen   '; 'Livingroom'; 'Office    '];
	rooms_desc = ['Sessions '];
    scene = struct('act', {});
    dtset = struct('tr', {}, 'te', {}, 'lb_tr', {}, 'lb_te', {});
    dts = struct('tr', {}, 'te', {}, 'lb_tr', {}, 'lb_te', {});
    featC = struct('activity', {}); %for combined features from sk1 and sk2
    test_featC = struct('activity', {}); 
    %1-handshake 2-hug 3-help walk 4-help stand-up 5-fight/punch 6-push 7-talk  8-draw attention
    
    % defining the proper activities per scenario according to CAD60
    scene(1).act = {1 2 3 4 5 6 7 8};
%    scene(2).act = {2 4 12 1 14};
%    scene(3).act = {4 10 11 12 1 14};
%    scene(4).act = {2 4 8 9 1 14};
%    scene(5).act = {2 3 4 13 1 14};
    
    tic;
    fprintf('\n Extracting Features...');
 
    [feat, feat2] = sk_features2(sk_matfile);
    [test_feat, test_feat2] = mono_sk_features2(sk_test_matfile);
    save('/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/Useful_files/individual_test_feat.mat','test_feat', 'test_feat2')
    save('/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/Useful_files/individual_tr_feat.mat','feat', 'feat2')
	load('/media/ccoppola/52FFC5D351F60A1C/Social-Activities-Code/Useful_files/k2_unsegmented.mat')
% 	% get features from individual skeletons 1 and 2 and combine them into a matrix for classification
	for l=1:num_persons
	   for c=1:numact
	       %maxtime=min(size(feat(l).activity{c},1),size(feat2(l).activity{c},1));
	       featC(l).activity{c} = [feat(l).activity{c}; feat2(l).activity{c}];
	       %feat(l).activity{c} = [feat(l).activity{c}(1:maxtime,:),feat2(l).activity{c}(1:maxtime,:)];
	   end
    end
    
    for l=1:size(test_feat,2)
	   for c=1:size(test_feat(1).activity,2)
	       %maxtime=min(size(feat(l).activity{c},1),size(feat2(l).activity{c},1));
	       test_featC(l).activity{c} = [test_feat(l).activity{c}; test_feat2(l).activity{c}];
	       %feat(l).activity{c} = [feat(l).activity{c}(1:maxtime,:),feat2(l).activity{c}(1:maxtime,:)];
	   end
    end
% 	
%     save('social_featSk1.mat', 'feat','-v7.3');
%     save('social_featSk2.mat', 'feat2','-v7.3');
%     save('social_featSk1_2.mat', 'featC','-v7.3');

% load('social_featSk1_2.mat');
% load('social_featSk2.mat');

    fprintf('  Done.\n');
    toc
    
    format short;
    
 
    fprintf('\nPreparing Data (spliting training and test sets)...\n');
 
    for r=1:size(rooms_desc, 1)
        
         fprintf('\n %s\n', rooms_desc(r, :));
         fprintf('===================\n');        
        
        
         % preparing the features of CAD60 in order to have each environment
         % and the respective activities - preapring the data for one-leave out
         % cross validation setup

         for i=1:size(test_featC,2)
%                % selection / inverting the learning and test => learn p2, p3, p4 and test p1; learn p3,p4,p1 and test p2; and so on... 
%                n = 1;
%                for k=1:num_persons
%                    if (k ~= i)
%                       t(n) = k; 
%                       n = n+1;
%                    end
%                end
% 
% 
%                fprintf('\nOrganizing Data for Learning Persons: ');
%                display(t);
%                fprintf('\nfor testing on Person: %d \n', i);
%                fprintf('------------------------\n\n');
% 
%                % getting each activity corresponding to each scenario
%                % performed by 3 persons as training set               
%                [Xtrain, Xtest] = getTrainTestData(r, scene, featC, t, i, test_featC);
%                


               Tr=(squeeze(struct2cell(featC)));
               Tr=cat(1,Tr{:});
               for k=1:numact
                    
                    XTrain{k}=cat(1,Tr{:,k});
                    labTrain{k}=k.*ones(size(cat(1,XTrain{k}),1),1);
                    
               end
                
               Te=(squeeze(struct2cell(test_feat)));
               XTest=cat(1,Te{:});
               labTest=Class_Truth;

               Te2=(squeeze(struct2cell(test_feat2)));
               XTest2=cat(1,Te2{:});
               
               %nomalization
               [ANorm_tr, ANorm_te] = apply_norm([cat(1,XTrain{:}),cat(1,labTrain{:})], [XTest{:},labTest]);
               [~, ANorm_te2] = apply_norm([cat(1,XTrain{:}),cat(1,labTrain{:})], [XTest2{:},labTest]);
               save('Normalized_Individual_Test_Features','ANorm_te','ANorm_te2');
% 			   dtset(r).lb_tr{i} = Xtrain(:, end);
% 			   dtset(r).lb_te{i} = Xtest(:, end);               
% 			   dtset(r).tr{i}  =  ANorm_tr; % Xtrain(:, 1:end-1); %   
% 			   dtset(r).te{i} = ANorm_te; % Xtest(:, 1:end-1);  % skeleton 1
%                
               
               dtset(r).tr{i} =    ANorm_tr; % Xtrain(:, 1:end-1); %   
			   dtset(r).te{i} =    ANorm_te; % Xtest(:, 1:end-1);  % 
			   dtset(r).lb_tr{i} = cat(1,labTrain{:});
			   dtset(r).lb_te{i} = labTest;
                        
               
               
%                %no normalized set
%                dts(r).tr{i} = Xtrain(:, 1:end-1);
%                dts(r).te{i} = Xtest(:, 1:end-1);
%                dts(r).lb_tr{i} = Xtrain(:, end);
%                dts(r).lb_te{i} = Xtest(:, end);                  
%             
%            
 
         end

         
    end
       
    
      
    disp('Done.');

    fprintf('\n\n Calling the Classification Stage ...\n');
    % calling base classifier 1-NB and 2-ANN, 4 (persons) tests per room  
    k2_classify(dtset, [1 2 3], num_tests);
    
    toc
end
        
    
    
    
    % Prepare the training set per scenario with the corresponding
    % activities, given the persons to learn
    % input: 
    %      room -> index corresponding the scenarios: 1-bath; 2-bed;
    %  3-kitchen; 4-living room; 5-office
    %      feat -> struct containing all features data from N persons and X
    %  scenarios with Z activities
    %      t -> a vector containing the persons index to learn, e.g. 2,3,4...
    %  (because the test will be the person 1) and so on...
    %      i -> an integer value indicating if the test person is 1, 2, 3...
    %
    % output:
    %      Xtrain-> a matrix of frames-by-N (columns 1:N=features; last column=labels)
    function [Xtrain, Xtest] = getTrainTestData(room, scene, featC, t, i, feat)
              
       % getting each activity corresponding to each scenario
       % performed by 3 persons as training set

        Xtrain = [];
        Xtest  = [];
        for x=1:size(scene(room).act,2)
            
          for y=1:size(t,2)
              lb = [];
              lb(1:size(featC(t(y)).activity{scene(room).act{x}},1),1) = scene(room).act{x};
              Xtrain = [Xtrain; featC(t(y)).activity{scene(room).act{x}} lb];
          end

         %test set
         lb = [];
         lb(1:size(feat(i).activity{scene(room).act{x}},1),1) = scene(room).act{x};
         Xtest = [Xtest; feat(i).activity{scene(room).act{x}} lb];    

        end                  
       
    end
        
        
     % Function to apply a normalization in the training and test set (per class / labels)
     %   input: Mtr (training) a matrix of nxm where the n lines are the frames and the
     %          column m are the features types. The last column is the
     %          label.
     %          Mte (test) similar to Mtr
     %
     %   output: ANorm_tr and ANorm_te the normalized feature matrices
     function [ANorm_tr, ANorm_te] = apply_norm(Mtr, Mte)
         
         ANorm_tr = [];
         ANorm_te = [];
         A = [];
         B = [];
         
         l = unique(Mtr(:,end));
         maxclass = [];
         minclass=[];
         for i=1:size(l,1)
             A=[];
             A = Mtr(Mtr(:,end) == l(i), 1:end-1);
             maxclass(i, :) = max(A);
             minclass(i, :) = min(A);
         end
         
         maxval = max(maxclass);
         minval = min(minclass);
         ranges = maxval - minval;
         A=[];
         A = Mtr(:, 1:end-1);
         B=[];         
         B = Mte(:, 1:end-1);
      
         
         ANorm_tr  = (A - repmat(minval, size(A,1), 1)) ./ repmat(ranges, size(A,1), 1);
         ANorm_te  = (B - repmat(minval, size(B,1), 1)) ./ repmat(ranges, size(B,1), 1);

     end 
         
         
   function [ftrain, ftest] = apply_filt(Xtrain, Xtest)
     for k=1:(size(Xtrain,1))
        ftrain(k,:) = medfilt1(Xtrain(k,:), 4);
     end
     for k=1:(size(Xtest,1))
        ftest(k,:) = medfilt1(Xtest(k,:), 4);
     end     
   end
             
         
         
