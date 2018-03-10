%[ skeleton_bags ] = k2_skeleton_text_to_matlab( '/home/ccoppola/Social_Activity_Detection_Data/Output' );
%[ annotation] = k2_csv_annotationtostruct('/home/ccoppola/Social_Activity_Detection_Data/Annotation.csv');
[ skeleton_bags ] = k2_skeleton_text_to_matlab( datafolder );
[ annotation] = k2_csv_annotationtostruct(annotationpath);


G_Truth=annotation;
for k=1:length(annotation)
    
    G_Truth(k).segmentation_user1=any(skeleton_bags{k,2}>annotation(k).segmentation_user1(1,:)&skeleton_bags{k,2}<annotation(k).segmentation_user1(2,:),2);
    G_Truth(k).segmentation_user2=any(skeleton_bags{k,2}>annotation(k).segmentation_user2(1,:)&skeleton_bags{k,2}<annotation(k).segmentation_user2(2,:),2);
    
end



clear k annotation;