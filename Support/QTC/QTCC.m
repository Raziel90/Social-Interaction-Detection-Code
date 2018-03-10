function [q1 q2 q3 q4 qd] = QTCC(x1, y1, x2, y2)
% QTC double-Cross
% x1..y2 are the vectors of coordinates of points 1 and 2
% q1..q4 are the vectors of QTC symbols (q1, q2: distance relations; q3, q4 side relations)
% qd is the vetor of qualitative diatnces between adjacent QTC states
% IMPORTANT: the QTC symbols for the first and and last coordinates are undefined, therefore q1..q4 refer to the second set of coordinates, in the x1..y2 vectors, until the second-last.
len = length(x1);
if (len < 3 | len ~= length(y1) | len ~= length(x2) | len ~= length(y2))
    error('Input vectors must have same lenght >= 3');
end

% distance relations
d1a = pdist2(x1(1:end-2), y1(1:end-2), x2(2:end-1), y2(2:end-1));
d1b = pdist2(x1(2:end-1), y1(2:end-1), x2(2:end-1), y2(2:end-1));
d1c = pdist2(x1(3:end), y1(3:end), x2(2:end-1), y2(2:end-1));

q1 = (d1a < d1b & d1b < d1c) + -1 * (d1a > d1b & d1b > d1c);

d2a = pdist2(x2(1:end-2), y2(1:end-2), x1(2:end-1), y1(2:end-1));
d2b = pdist2(x2(2:end-1), y2(2:end-1), x1(2:end-1), y1(2:end-1));
d2c = pdist2(x2(3:end), y2(3:end), x1(2:end-1), y1(2:end-1));

q2 = (d2a < d2b & d2b < d2c) + -1 * (d2a > d2b & d2b > d2c);

%side relations
r1a = [x1(1:end-2); y1(1:end-2)] - [x1(2:end-1); y1(2:end-1)];
r1b = [x2(2:end-1); y2(2:end-1)] - [x1(2:end-1); y1(2:end-1)];
r1c = [x1(3:end); y1(3:end)] - [x1(2:end-1); y1(2:end-1)];
s1a = side(r1a, r1b);
s1c = side(r1c, r1b);

q3 = (s1a < 0 & s1c > 0) + -1 * (s1a > 0 & s1c < 0);
 
r2a = [x2(1:end-2); y2(1:end-2)] - [x2(2:end-1); y2(2:end-1)];
r2b = [x1(2:end-1); y1(2:end-1)] - [x2(2:end-1); y2(2:end-1)];
r2c = [x2(3:end); y2(3:end)] - [x2(2:end-1); y2(2:end-1)];
s2a = side(r2a, r2b);
s2c = side(r2c, r2b);

q4 = (s2a < 0 & s2c > 0) + -1 * (s2a > 0 & s2c < 0);

% qualitative distance between adjacent QTC states
qd = abs(q1(2:end)-q1(1:end-1)) + abs(q2(2:end)-q2(1:end-1)) + abs(q3(2:end)-q3(1:end-1)) + abs(q4(2:end)-q4(1:end-1));

end


function d2 = pdist2(x1, y1, x2, y2)
d2 = (x1-x2).^2 + (y1-y2).^2;
end

function s = side(a, b)
cp= cross([a; zeros(1,length(a))], [b; zeros(1,length(b))]);
s = cp(3,:);
end
