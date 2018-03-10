function [  ] = k2_skeleton_view( skeleton_sample,color_string,marker_size,Line_width,labels_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%skelv2=reshape(processed_skeletons{1,3}(1).data(50,:)',7,25)';
%bone_pos=[0,1;1,20;20,2;2,3;20,4;20,8;8,9;9,10;10,11;11,23;10,24;4,5;5,6;6,7;7,21;6,22;0,16;0,12;16,17;17,18;18,19;12,13;13,14;14,15]+1;
bone_pos=[5,4;4,3;3,2;2,1;1,14;1,18;14,15;15,16;16,17;18,19;19,20;20,21;3,6;6,7;7,8;8,9;9,22;8,23;3,10;10,11;11,12;12,13;13,24;12,25];
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
dx=0.03;
dy=0.04;
c=cellstr(num2str((1:25)'))'; 
skelv2=reshape(skeleton_sample,7,25)';
figure;
hold on
for b=1:size(bone_pos,1)
   g=line(skelv2(bone_pos(b,:),1),skelv2(bone_pos(b,:),2),skelv2(bone_pos(b,:),3),'LineWidth',Line_width);
   uistack(g,'bottom')
end

scatter(skelv2(:,1),skelv2(:,2),marker_size*ones(size(skelv2(:,2))),['.',color_string]);

if labels_flag==1
    textpos1=skelv2(:,1)+dx;
    textpos1([6:9])=textpos1([6:9])-3*dx;
    textpos2=skelv2(:,2)+dy;
    text(textpos1, textpos2, c)
end


hold off
xlim([2.5 3.1])
ylim([-1.5 1])

axis equal
% axis([2.5 3.1 -1.5 1])
axis off
end

