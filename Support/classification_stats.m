function [ confmat,acc_prec_rec,fig_h] = classification_stats( labels,estimated_labels,plotit,normalizeit,label_names)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
labels=reshape(labels,1,length(labels));
estimated_labels=reshape(estimated_labels,1,length(estimated_labels));
class_set=unique([labels(:);estimated_labels(:)]');
slack=max([sum(class_set==0),0]);

stats=zeros(length(class_set),4);
confmat=zeros(length(class_set));
precision=zeros(length(class_set),1);
recall=zeros(length(class_set),1);

for c=1:length(class_set)
    
   stats(c,1)=sum(estimated_labels==(c-slack)&labels==(c-slack)); %TP
   stats(c,2)=sum(estimated_labels==(c-slack)&labels~=(c-slack)); %FP
   stats(c,3)=sum(estimated_labels~=(c-slack)&labels~=(c-slack)); %TN
   stats(c,4)=sum(estimated_labels~=(c-slack)&labels==(c-slack)); %FN
%    confmat(c,c)
accuracy(c)=sum(stats(c,[1,3]))/sum(stats(c,:));
precision(c)=stats(c,1)/sum(stats(c,[1,2]));
recall(c)=stats(c,1)/sum(stats(c,[1,4]));
end

%     [c,confmat]=confusion(estimated_labels,labels);
    confmat = confusionmat(labels,estimated_labels);
    acc_prec_rec=[accuracy',precision,recall];
    normconfmat=confmat./sum(confmat,2).*100;
    if nargin>3&&normalizeit==1
        confmat=normconfmat;
        formatstr='%0.2f';
    else
        formatstr='%0.0f';
    end
    if nargin>2 && plotit==1
    h=imagesc(normconfmat);            %# Create a colored plot of the matrix values
    colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                             %#   black and lower values are white)

    textStrings = num2str(confmat(:),formatstr);  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:length(class_set));   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(imgca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(normconfmat(:) > midValue,1,3);  %# Choose white or black for the
                                                 %#   text color of the strings so
                                                 %#   they can be easily seen over
                                                 %#   the background color
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    if nargin > 4
    set(imgca,'XTick',1:length(label_names),...                         %# Change the axes tick marks
            'XTickLabel',label_names,...  %#   and tick labels
            'YTick',1:length(label_names),...
            'YTickLabel',label_names,...
            'YTickLabelRotation',0,...
            'TickLength',[0 0],...
            'XTickLabelRotation',-45)
            %'XAxisLocation','Top'...
            %);
	


    end
end
if nargout > 2
 fig_h=imgcf;
else
    fig_h=[];
end






end

