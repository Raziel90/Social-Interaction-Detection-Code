function [qtcidx] = qtc2idx(qtcc,strlen)
if nargin<2
    strlen=4;
end
positional_values=3.^(strlen-1:-1:0);
%iqtc=qtcc+1;
qtcidx=(sum(repmat(positional_values,size(qtcc,1),1).*(qtcc+1),2)+1);
end