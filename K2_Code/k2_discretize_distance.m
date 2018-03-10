function [distance_states] = k2_discretize_distance(skeleton_bags,min_distances_sector,percentage_of_slack,joint_idx)
% proxemic_sectors=[0.45,1.2,3.7,7.6];
% min_distances_sector=[0.15,0.46,0.76,1.22,2.1,3.7,7.6];
% qtc_sequences=cell(size(skeleton_bags,1),1);
numdimensions=7;
numjoints=25;
if length(min_distances_sector)>2
incoming_th=min_distances_sector-[0,diff(min_distances_sector)]*percentage_of_slack(1);
outgoing_th=min_distances_sector+[0,diff(min_distances_sector)]*percentage_of_slack(2);

for k=1:size(skeleton_bags,1)
    
    lenclip=length(skeleton_bags{k,2});
    skel1=nan(numdimensions,numjoints,lenclip);
    skel2=nan(numdimensions,numjoints,lenclip);
    skel1(:,:,skeleton_bags{k,3}(1).temporal_presence)=reshape(skeleton_bags{k,3}(1).data',numdimensions,numjoints,length(skeleton_bags{k,3}(1).data));
    skel2(:,:,skeleton_bags{k,3}(2).temporal_presence)=reshape(skeleton_bags{k,3}(2).data',numdimensions,numjoints,length(skeleton_bags{k,3}(2).data));
    skel1=skel1([1,3],joint_idx,:);
    skel2=skel2([1,3],joint_idx,:);
%     distances=sqrt(sum(squeeze(skel1-skel2).^2));
%     distances=smooth(distances,20,'sgolay');
    distances=medfilt1(sqrt(sum(squeeze(skel1-skel2).^2)),20);
    distance_states{k}=zeros(size(distances));
%     distance_states{n}(distance<=min_distances_sector(1))=0
    for n=1:(length(min_distances_sector)-1)
%         distance_states{k}(distances>min_distances_sector(n)&distances<=min_distances_sector(n+1))=n;
        distance_states{k}(distances>incoming_th(n)&distances<=outgoing_th(n+1))=n;
    end
    distance_states{k}(distances>min_distances_sector(n+1))=n+1;
        
end
else
    distance_states=[];
end

end