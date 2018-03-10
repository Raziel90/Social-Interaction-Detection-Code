function [  ] = NITE2_skeleton_view( skeleton_sample,color_string,marker_size,Line_width,labels_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%skelv2=reshape(processed_skeletons{1,3}(1).data(50,:)',7,25)';
bone_pos=[1,2;2,4;6,2;4,3;3,6;3,8;3,10;8,10;10,11;11,15;8,9;9,14;4,5;5,12;6,7;7,13];
if nargin<2
   color_string=''; 
end
if nargin<3
   marker_size=500; 
end
if nargin<4
   Line_width=1; 
end
if nargin<5
   labels_flag=1; 
end
dx=0.04;
dy=0.02;
c=cellstr(num2str((1:15)'))'; 
skelv2=reshape(skeleton_sample',6,15)';
figure;
hold on
for b=1:size(bone_pos,1)
   g=line(skelv2(bone_pos(b,:),1),skelv2(bone_pos(b,:),3),skelv2(bone_pos(b,:),2),'LineWidth',Line_width);
   uistack(g,'bottom')
end
scatter(skelv2(:,1),skelv2(:,3),marker_size*ones(size(skelv2(:,2))),['.',color_string])
if labels_flag==1
    textpos1=skelv2(:,1)+dx;
    textpos1([4,5,12])=textpos1([4,5,12])-3.5*dx;
    textpos2=skelv2(:,3)+dy;
    text(textpos1,textpos2, c)
end
% xlim([2.5 3.1])
% ylim([-1.5 1])

axis equal
% axis([2.5 3.1 -1.5 1])
axis off
end

