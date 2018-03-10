function [ output_args ] = QTC_numseq2str( QTC_seq )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
cellfun(@(Y)strrep(strrep(strrep(Y','-1','-'),'1','+'),' ','')',mat2cell(num2str(input_sequences{1}')',size(num2str(input_sequences{1}')',1),ones(1,size(input_sequences{1},2))),'UniformOutput',false)

cell2mat(cellfun(@(Y)strrep(strrep(strrep(Y','-1','-'),'1','+'),' ','')',mat2cell(num2str(input_sequences{1}')',size(num2str(input_sequences{1}')',1),ones(1,size(input_sequences{1},2))),'UniformOutput',false))'

% QTC_NUM_STR=@(QTC_seq)cellfun(@(X)strrep(strrep(strrep(num2str(X),'-1','-'),'1','+'),' ',''),mat2cell(QTC_seq,ones(size(QTC_seq,1),1),4),'UniformOutput',false);
QTC_NUM_STR=@(QTC_seq)cellfun(@(X)strrep(strrep(strrep(num2str(X),'-1','-'),'1','+'),' ',''),QTC_seq,'UniformOutput',false);
end

