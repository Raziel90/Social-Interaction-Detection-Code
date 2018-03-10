%================================================================================
% ISR-UC - Institute of Systems and Robotics, University of Coimbra, Portugal
%================================================================================
% Authors: Diego R. Faria, Cristiano Premebida  -  Last update: November, 2014
%
% Description: Naive Bayes with Gaussian distribution - base classifier for the DBMM - generative classifier
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% input:          X_train -> frames x N vector, training set with the features for all classes, sequentially, e.g. 1:100 class1; 101:300 class2; 300:N class3 ... 
%                  X_test -> similar to X_train but with the test set features (it can be just one single class for testing)
%                 label_tr-> a vector (framesx1) with the labels for each class in the training set
%                 label_te-> similar to label_tr, but using the test set labels
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
% output: classification results - Rtr -> a matrix (with dimension Frames x Num of classes) containing the classif. results (probabilities) for the training set
%                                  Rte -> similar to Rtr, but with results for the test set ...
%                                  cm_tr_nb -> confusion matrix for the training set (square matrix)
%                                  cm_te_nb -> confusion matrix for the test set (square matrix)
%                                  learned  -> learned likelihood (training parameters)
%------------------------------------------------------------------------------------------------------------------------------------------------------------


function [Ytr, Yte, cm_tr_svm, cm_te_svm, learned] = SVM_tr_te(X_train, X_test, label_tr, label_te, ker, cost, gamma)
   
    % ensure that the label vectors are with dimensions frames x 1
    if (size(label_tr, 2) > size(label_tr, 1))
       label_tr = label_tr';
    end
    if (size(label_te, 2) > size(label_te, 1))
       label_te = label_te';
    end

	size_class = size(unique(label_tr),1);

	%get the boundaries (transitions) between classes in the training and test set
	label_bound_tr = [];
    classes_tr = unique(label_tr);
	for c=1:size_class
		ini =  find(label_tr(:,1)==classes_tr(c), 1, 'first');
		fim =  find(label_tr(:,1)==classes_tr(c), 1, 'last');
		label_bound_tr = [label_bound_tr; ini fim];
    end 
    
	label_bound_te = [];
    classes_te = unique(label_te);
	for c=1:size_class
		ini =  find(label_te(:,1)==classes_te(c), 1, 'first');
		fim =  find(label_te(:,1)==classes_te(c), 1, 'last');
		label_bound_te = [label_bound_te; ini fim];
    end 
    
    if nargin<7
       gamma=1/min(size(X_train));
       if nargin<6
           cost=1;
           if nargin<5
               ker=0;
           end
       end
       
    end

    acc_prec_recall_svm =[];

    learned = libsvmtrain(label_tr, X_train, ['-s 0 -t ',num2str(ker),' -c ',num2str(cost),' -g ',num2str(gamma),' -e 0.01 -m 4000 -b 1 -q']);    % train-c 1000 -g 10 -b 1
    
    %test on the training set
   
    [predicted_label, accuracy, Ytr] = libsvmpredict(label_tr, X_train, learned, '-b 1');
    [predicted_label, accuracy, Yte] = libsvmpredict(label_te, X_test, learned, '-b 1');   % faz a predicao    
    
    

	% get the confusion matrix
	cm_te_svm = count_correct(Yte', label_te);
	cm_tr_svm = count_correct(Ytr', label_tr);
	cm_te_svm = cm_te_svm';
	cm_tr_svm = cm_tr_svm';
    
    
	 % class_label = [];
	 % for c=1:size_class
		 % class_label = [class_label; label_bound_te(c,1) label_bound_te(c,2)];
	 % end
	 % name = [dtset_name '_label.mat'];
	 % save(name, 'class_label');
	 
	 % class_labeltr = [];
	 % for c=1:size_class
		 % class_labeltr = [class_labeltr; label_bound_tr(c,1) label_bound_tr(c,2)];
	 % end
	 % name = [dtset_name '_train_label.mat'];
	 % save(name, 'class_labeltr');

end




 % count the correct classification (accuracy) after the frame by frame classification
 % input: pclass = a matrix of dimension number of classifiers (line) by number of frames (column), containing the classification probabilities
 %         class_label =  a matrix of Cx2 (number of labels (classes) by column 1 initial and column 2 final position of the labels (nº frames))
 %                        example: clabel = [1 100; 101 200; 201 500] - each line correspond to the number of classes, in this case 3, 
 %                        and columns the initial and final position of each label regarding the frame numbers
 % output: corr = confusion matrix with dimension classes x classes
function corr = count_correct(pclass, class_label)
   corr=[];
   % count the error for the true class
   present_classes=sort(unique(class_label));
   num_class=length(unique(class_label));
   for c=1:num_class
	   classif = zeros(1,num_class);
       class_idx=find(present_classes(c)==class_label);
	   for m=1:length(class_idx)
		   [~,indi] = sort(pclass(:,class_idx(m)),'descend');
		   classif(indi(1)) = classif(indi(1)) + 1;            
	   end
	   classif=classif./length(class_idx)*100;
	   corr(:,c) = classif';
   end
end
