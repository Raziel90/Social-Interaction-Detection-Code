function [a12,a21] = k2_compute_body_orientation(skeletons,J1,J2,plot )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin <4
    plot=false;
end

numdimensions=7;
numjoints=25;
a12=cell(size(skeletons,1),1);
a21=cell(size(skeletons,1),1);
for k=1:length(skeletons)
    lenclip=length(skeletons{k,2});
    skel1=nan(numdimensions,numjoints,lenclip);
    skel2=nan(numdimensions,numjoints,lenclip);
    skel1(:,:,skeletons{k,3}(1).temporal_presence)=reshape(skeletons{k,3}(1).data',numdimensions,numjoints,length(skeletons{k,3}(1).data));
    skel2(:,:,skeletons{k,3}(2).temporal_presence)=reshape(skeletons{k,3}(2).data',numdimensions,numjoints,length(skeletons{k,3}(2).data));
    
    t1=squeeze(skel1([1,2,3],[J1],:)-skel1([1,2,3],[J2],:))';
    t2=squeeze(skel2([1,2,3],[J1],:)-skel2([1,2,3],[J2],:))';
    
    zt1=find(t1(:,1)==0);
    zt2=find(t2(:,1)==0);
    
    for n=1:length(zt1)
        t1(zt1(n),:)=t1(zt1(n)-1,:);
    end
    
    for n=1:length(zt2)
        t2(zt2(n),:)=t2(zt2(n)-1,:);
    end
    
    
    
    n1=[-t1(:,3),-t1(:,1)];
    n2=[-t2(:,3),t2(:,1)];
    n1=n1./sqrt(sum(n1.^2,2));
    n2=n2./sqrt(sum(n2.^2,2));
    
    m1=squeeze(mean(skel1([1,2,3],[J1,J2],:),2))';
    m2=squeeze(mean(skel2([1,2,3],[J1,J2],:),2))';
    
    v12=m2(:,[1,3])-m1(:,[1,3]);
    v12=v12./sqrt(sum(v12.^2,2));
    v21=-v12;
    a12{k}=acos(dot(n1,v12,2))';
    a21{k}=acos(dot(n2,v21,2))';
    if plot == true
        for h=1:lenclip
%             plot3(squeeze(skel1([1],[J1,J2],h)),squeeze(skel1([2],[J1,J2],h)),squeeze(skel1([3],[J1,J2],h))),hold on,scatter([squeeze(skel1([1],[J1,J2],h)),m1(h,1)],[squeeze(skel1([2],[J1,J2],h)),m1(h,2)]),xlim([-1,1]),quiver3(m1(h,1),m1(h,2),m1(h,3),n1(h,1),0,n2(h,2))
%             plot3(squeeze(skel2([1],[J1,J2],h)),squeeze(skel2([2],[J1,J2],h)),squeeze(skel2([3],[J1,J2],h))),hold on,scatter([squeeze(skel2([1],[J1,J2],h)),m2(h,1)],[squeeze(skel2([2],[J1,J2],h)),m2(h,2)]),quiver3(m2(h,1),m2(h,2),m2(h,3),n2(h,1),0,n2(h,2))
%             quiver(m1(h,1),m1(h,3),n1(h,1),n1(h,2)),hold on,quiver(m2(h,1),m2(h,3),n2(h,1),n2(h,2))%,ylim([-1,1]),xlim([-1,1])
            quiver([m1(h,1);m1(h,1)],[m1(h,3);m1(h,3)],[n1(h,1);v12(h,1)],[n1(h,2);v12(h,2)]),hold on,quiver([m2(h,1);m2(h,1)],[m2(h,3);m2(h,3)],[n2(h,1);v21(h,1)],[n2(h,2);v21(h,2)])
    %         
    %         plot3(squeeze(skel1([1],[J1,J2],h)),squeeze(skel1([2],[J1,J2],h)),squeeze(skel1([3],[J1,J2],h))),hold on,scatter3([squeeze(skel1([1],:,h)),m1(h,1)],[squeeze(skel1([2],:,h)),m1(h,2)],[squeeze(skel1([3],:,h)),m1(h,3)]),xlim([-1,1]),quiver3(m1(h,1),m1(h,2),m1(h,3),n1(h,1),0,n2(h,2))
    %         plot3(squeeze(skel2([1],[J1,J2],h)),squeeze(skel2([2],[J1,J2],h)),squeeze(skel2([3],[J1,J2],h))),hold on,scatter3([squeeze(skel2([1],:,h)),m2(h,1)],[squeeze(skel2([2],:,h)),m2(h,2)],[squeeze(skel2([3],:,h)),m2(h,3)]),ylim([-1,1]),quiver3(m2(h,1),m2(h,2),m2(h,3),n2(h,1),0,n2(h,2))
    %         view([0,88])

            pause(1/15)
            clf
        end
    end
end

