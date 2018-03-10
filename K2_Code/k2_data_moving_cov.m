function [ data_moving_cov ] = k2_data_moving_cov( input_data,window,dimensions,eig_mode)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
num_dim=length(dimensions);
if nargin<4
    eig_mode=1;
    if nargin<3
        dimensions=[2,3];
        if nargin<2
            window=15;
        end
    end
    
end


data_moving_cov=cell(size(input_data));
for h=1:length(input_data);
        
    if eig_mode==0
        
        data_moving_cov{h}=zeros(num_dim.^2,max(size(input_data{h})));
    else
        data_moving_cov{h}=zeros(1,max(size(input_data{h})));
    end

m_cov=zeros(num_dim,num_dim,max(size(input_data{h}(1,:))));

for j=2:(window-1)
m_cov(:,:,j)=cov(input_data{h}(dimensions,1:j)');    
end
m_cov(:,:,window:end)=moving_cov(input_data{h}(dimensions,:)',window);

for n=2:size(m_cov,3)
    
    if eig_mode==0
        
        data_moving_cov{h}(:,n)=reshape(logm(m_cov(:,:,n)+1),1,num_dim.^2);
    else
        data_moving_cov{h}(:,n)=max(eig(m_cov(:,:,n)));
    end
    
    
end


end





end

