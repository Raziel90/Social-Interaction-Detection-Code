function [ annotation] = k2_csv_annotationtostruct( annotationpath )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
annotationpath='/home/ccoppola/Social_Activity_Detection_Data/Annotation.csv';
end


[Annotation1] = importfile(annotationpath);



numfiles=sum(~cell2mat(cellfun(@(X)contains(string(X),'user_'),Annotation1(:,2),'UniformOutput',false)));
Annotation1(:,size(Annotation1,2))=[];
filenamepositions=find(~cell2mat(cellfun(@(X)contains(string(X),'user_'),Annotation1(:,2),'UniformOutput',false)));
filenames=string(Annotation1(filenamepositions,2));
sessionnames=reshape(repmat(string(Annotation1(filenamepositions(1:2:length(filenamepositions)),1)),1,2)',1,numfiles)';
annotation=repmat(struct('session_name','','namefile','','segmentation_user1',[],'segmentation_user2',[]),20,1);

for k=1:numfiles
    annotation(k).session_name=sessionnames(k);
    annotation(k).namefile=filenames(k);
    nonnan=cellfun(@(X)~isnan(str2double(X)),Annotation1(filenamepositions(k)+1,1:end));
    annotation(k).segmentation_user1=cellfun(@str2double,Annotation1(filenamepositions(k)+1,nonnan))';
    annotation(k).segmentation_user1=reshape(annotation(k).segmentation_user1,2,max(size(annotation(k).segmentation_user1))/2);
    annotation(k).segmentation_user2=cellfun(@str2double,Annotation1(filenamepositions(k)+2,nonnan))';
    annotation(k).segmentation_user2=reshape(annotation(k).segmentation_user2,2,max(size(annotation(k).segmentation_user2))/2);
    
end



end









function [output] = importfile(filename)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [SESSION_1,RECORDING_1482159915,START,STOP,START1,STOP1,START2,STOP2,START3,STOP3]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [SESSION_1,RECORDING_1482159915,START,STOP,START1,STOP1,START2,STOP2,START3,STOP3]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [session_1,recording_1482159915,start,stop,start1,stop1,start2,stop2,start3,stop3] = importfile('Annotation.csv',2, 30);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/01/17 16:37:52

%% Initialize variables.
delimiter = ',';

startRow = 1;
endRow = inf;


%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
output = horzcat(dataArray{:,:});
% stop = dataArray{:, 4};
% start1 = dataArray{:, 5};
% stop1 = dataArray{:, 6};
% start2 = dataArray{:, 7};
% stop2 = dataArray{:, 8};
% start3 = dataArray{:, 9};
% stop3 = dataArray{:, 10};
end