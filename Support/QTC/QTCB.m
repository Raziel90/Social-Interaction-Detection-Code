function [q1 q2 qd] = QTCB(x1, y1, x2, y2)
% QTC Basic
% Input
%  x1, y1 : cordinates of point 1
%  x2, y2 : cordinates of point 2
% Output
%  q1 : QTC_B symbol of point 1 wrt point 2
%  q2 : QTC_B symbol of point 2 wrt point 1
%  qd : qualitative distance between adjacent QTC_B states

len = length(x1);
if (len < 3 | len ~= length(y1) | len ~= length(x2) | len ~= length(y2))
    error('Inout vectors must have same lenght >= 3');
end

d1a = pdist2(x1(1:end-2), y1(1:end-2), x2(2:end-1), y2(2:end-1));
d1b = pdist2(x1(2:end-1), y1(2:end-1), x2(2:end-1), y2(2:end-1));
d1c = pdist2(x1(3:end), y1(3:end), x2(2:end-1), y2(2:end-1));

q1 = (d1a < d1b & d1b < d1c) + -1 * (d1a > d1b & d1b > d1c);

d2a = pdist2(x2(1:end-2), y2(1:end-2), x1(2:end-1), y1(2:end-1));
d2b = pdist2(x2(2:end-1), y2(2:end-1), x1(2:end-1), y1(2:end-1));
d2c = pdist2(x2(3:end), y2(3:end), x1(2:end-1), y1(2:end-1));

q2 = (d2a < d2b & d2b < d2c) + -1 * (d2a > d2b & d2b > d2c);

% distance between adjacent QTC states
qd = abs(q1(2:end)-q1(1:end-1)) + abs(q2(2:end)-q2(1:end-1));

end


function d2 = pdist2(x1, y1, x2, y2)
d2 = (x1-x2).^2 + (y1-y2).^2;
end
