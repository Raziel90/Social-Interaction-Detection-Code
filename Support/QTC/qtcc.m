function [ qtcc_rep, frames ] = qtcc( human, robot, accuracy,repetitions, plot )
%QTCC Summary of this function goes here
%   Creates (and plots) a QTC_C sequence from a given vector of points.
%   human and robot are x,y column vectors, accuracy is the distance to the
%   lines at which a point is considered to lay on the line.
%   plot = 1 turns on visualisation and saves them as framjes in frame to
%   be replayed afterwards.

qtcc_rep_tmp = [0 0 0 0];

end_idx = size(human,1);
if size(human,1) ~= size(robot,1)
    end_idx = min([size(human,1) size(robot,1)]);
end

for i=2:end_idx
    RL = [human(i-1,1:2); robot(i-1,1:2)];
    tmp1 = qtcTranslate(human(i-1,1:2), .5*(human(i-1,1:2)-robot(i-1,1:2)));
    tmp2 = qtcTranslate(robot(i-1,1:2), .5*(robot(i-1,1:2)-human(i-1,1:2)));
    RL_ext = [tmp1(1:2); tmp2(1:2)];
    rot_RL = orthogonalLine([human(i-1,1:2) robot(i-1,1:2)-human(i-1,1:2)], human(i-1,1:2));
    trans_RL_h = qtcTranslate([rot_RL(1,1:2); rot_RL(1,3:4)], (rot_RL(1,1:2)-rot_RL(1,3:4))/2);
    trans_RL_r = qtcTranslate(trans_RL_h, (robot(i-1,1:2)-human(i-1,1:2)));
    
    % robot movement: k
    k = [Inf, Inf];
    
    x0 = robot(i,1);
    y0 = robot(i,2);
    x1 = trans_RL_r(2,1);
    y1 = trans_RL_r(2,2);
    x2 = trans_RL_r(1,1);
    y2 = trans_RL_r(1,2);
    test=(x0 - x1) * (y2 - y1) - (x2 - x1) * (y0 - y1);
    d = abs(det([trans_RL_r(2,1:2)-trans_RL_r(1,1:2);robot(i,1:2)-trans_RL_r(1,1:2)]))/norm(trans_RL_r(2,1:2)-trans_RL_r(1,1:2));
    if test > 0 & abs(d) > accuracy
        k(1,1) = -1;
    elseif test < 0 & abs(d) > accuracy
        k(1,1) = 1;
    else
        k(1,1) = 0;
    end

    x0 = robot(i,1);
    y0 = robot(i,2);
    x1 = RL_ext(2,1);
    y1 = RL_ext(2,2);
    x2 = RL_ext(1,1);
    y2 = RL_ext(1,2);
    test=(x0 - x1) * (y2 - y1) - (x2 - x1) * (y0 - y1);
    d = abs(det([RL_ext(2,:)-RL_ext(1,:);robot(i,1:2)-RL_ext(1,:)]))/norm(RL_ext(2,:)-RL_ext(1,:)); 
    if test < 0 & abs(d) > accuracy
        k(1,2) = -1;
    elseif test > 0 & abs(d) > accuracy
        k(1,2) = 1;
    else
        k(1,2) = 0;
    end
    
    % movement of human: l
    l = [Inf, Inf];

    x0 = human(i,1);
    y0 = human(i,2);
    x1 = trans_RL_h(1,1);
    y1 = trans_RL_h(1,2);
    x2 = trans_RL_h(2,1);
    y2 = trans_RL_h(2,2);
    test = (x0 - x1) * (y2 - y1) - (x2 - x1) * (y0 - y1);
    d = abs(det([trans_RL_h(2,1:2)-trans_RL_h(1,1:2);human(i,1:2)-trans_RL_h(1,1:2)]))/norm(trans_RL_h(2,1:2)-trans_RL_h(1,1:2));
    if test > 0 & abs(d) > accuracy
        l(1,1) = -1;
    elseif test < 0 & abs(d) > accuracy
        l(1,1) = 1;
    else
        l(1,1) = 0;
    end

    x0 = human(i,1);
    y0 = human(i,2);
    x1 = RL_ext(1,1);
    y1 = RL_ext(1,2);
    x2 = RL_ext(2,1);
    y2 = RL_ext(2,2);
    test = (x0 - x1) * (y2 - y1) - (x2 - x1) * (y0 - y1);
    d = abs(det([RL_ext(2,:)-RL_ext(1,:);human(i,1:2)-RL_ext(1,:)]))/norm(RL_ext(2,:)-RL_ext(1,:));
    if test < 0 & abs(d) > accuracy
        l(1,2) = -1;
    elseif test > 0 & abs(d) > accuracy
        l(1,2) = 1;
    else
        l(1,2) = 0;
    end
    
    % print
    qtcc_rep_tmp = [qtcc_rep_tmp; [l(1,1) k(1,1) l(1,2) k(1,2)]];
    
    % plot
    if nargin == 4 & plot == 1
        qtcPlotPoses(human(i-1:i,:), robot(i-1:i,:), 2, 1, 1)
        axis equal
        line_rl = line(RL_ext(:,1), RL_ext(:,2), 'Color', 'm');
        rl_l = line(trans_RL_h(:,1)', trans_RL_h(:,2)', 'Color', 'm');
        rl_r = line(trans_RL_r(:,1)', trans_RL_r(:,2)', 'Color', 'm');
        title(['QTC_C state transition: (', num2str(qtcc_rep_tmp(end-1,:)), ') -> (',num2str(qtcc_rep_tmp(end,:)), ')']);
        frames(i-1) = getframe;
        pause(0.01);
        delete(line_rl);
        delete(rl_l);
        delete(rl_r);
    else
        frames = struct;
    end
end

if isempty(qtcc_rep_tmp)
    disp('!!!Empty qtcc_rep!!!')
    qtcc_rep = [0 0 0 0];
    frames = struct;
    return;
end

if size(qtcc_rep_tmp,1) == 0
    disp('!!!Just one data point!!!')
    qtcc_rep = [0 0 0 0];
    frames = struct;
    return;
end
if repetitions == 0
    qtcc_rep = qtcc_rep_tmp(1,:);
    j = 1;
    for i=2:size(qtcc_rep_tmp,1)
        if(~isequal(qtcc_rep(j,:), qtcc_rep_tmp(i,:)))
            qtcc_rep = [qtcc_rep; qtcc_rep_tmp(i,:)];
            j = j + 1;
        end
    end
else
    qtcc_rep = qtcc_rep_tmp;
end
 
end




function res = orthogonalLine(line, point)
%ORTHOGONALLINE Create a line orthogonal to another one.
%
%   PERP = orthogonalLine(LINE, POINT);
%   Returns the line orthogonal to the line LINE and going through the
%   point given by POINT. Directed angle from LINE to PERP is pi/2.
%   LINE is given as [x0 y0 dx dy] and POINT is [xp yp].
%
%   See also:
%   lines2d, parallelLine
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   19/02/2004 added control for multiple lines and/or points


N = max(size(point, 1), size(line, 1));

if size(point, 1)>1
    res = point;
else
    res = ones(N, 1)*point;
end

if size(line, 1)>1
    res(:,3) = -line(:,4);
    res(:,4) = line(:,3);
else
    res(:,3) = -ones(N,1)*line(4);
    res(:,4) = ones(N,1)*line(3);
end

res(:,3) = res(:,1)+res(:,3);
res(:,4) = res(:,2)+res(:,4);
end






function [ pose_trans ] = qtcTranslate( pose, trans_vec )
%TRANSLATE Summary of this function goes here
%   Translating every 2D point in pose by trans_vec


pose_trans = pose;

for i = 1:size(pose_trans,1)
    pose_trans(i,1:2)=pose_trans(i,1:2)+trans_vec;
end

end


function [ fig ] = qtcPlotPoses( human, robot, new_plot, plot_start, plot_end )
%PLOTPOSES Summary of this function goes here
%   Detailed explanation goes here

ps = 1;
pe = 1;

if nargin < 3 || new_plot == 1
    fg=figure;
    set(fg,'OuterPosition',[2638 326 562 505]);
elseif nargin >= 3 & new_plot == 2
    refresh
    set(gcf,'OuterPosition',[2638 326 562 505]);
end
if nargin >= 4 & plot_start == 0
    ps = 0; 
end
if nargin >= 5 & plot_end == 0
    pe = 0; 
end

hold on

tmp1 = plot(human(:,1), human(:,2), 'b-', 'LineWidth', 2);
tmp2 = plot(robot(:,1), robot(:,2), 'r-', 'LineWidth', 2);
% if nargin < 3 || new_plot == 1
    h = [tmp1; tmp2];
    set(h(1),'Displayname','Person1')
    set(h(2),'Displayname','Person2')
    legend(h,'Location','north')
% end

if ps == 1
    plot(human(1,1), human(1,2), 'rs')
    plot(robot(1,1), robot(1,2), 'rs')
end
if pe == 1
    plot(human(size(human,1),1), human(size(human,1),2), 'kd')
    plot(robot(size(robot,1),1), robot(size(robot,1),2), 'kd')
end

end
