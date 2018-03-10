
function [cumulative]=accumulate_joint_energy(Venerg,w,joint_ind)
    if nargin<3
        joint_ind=1:15;
    end
    for c=1:size(Venerg,1)
        cumulative{c}=sum(Venerg{c,1}([joint_ind],w+1:end))+sum(Venerg{c,2}([joint_ind],w+1:end));
    end
end


