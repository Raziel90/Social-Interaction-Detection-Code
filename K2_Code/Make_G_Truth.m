function [ class_vector ] = Make_G_Truth( Seq_user1,Seq_user2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

skel1=reshape(Seq_user1',6,15,length(Seq_user1));
skel2=reshape(Seq_user2',6,15,length(Seq_user2));
skel1=skel1(1:3,:,:);
skel2=skel2(1:3,:,:);
class_vector=zeros(length(Seq_user2),1);
for frame=1:length(Seq_user2)
h=figure(1); hold on
scatter3(skel1(1,:,frame),skel1(2,:,frame),skel1(3,:,frame))
scatter3(skel2(1,:,frame),skel2(2,:,frame),skel2(3,:,frame))
x=[];
while isempty(x)
 x = str2num(input(['what is the class now? t=',num2str(frame),'/',num2str(length(Seq_user2)),' : '],'s'));
end

class_vector(frame)=x;
clf(h)




end

