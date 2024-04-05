function [frame_length, thr, estimated_gain_dist,permitted] = ...
    assign_frame_length(estimated_gain_dist,gain)
global wmin
u = (estimated_gain_dist(1,2:end));
u1 = fliplr(cumsum(u,'reverse'));
ind = find(u1>=wmin,1,'first');
temp = u1(ind);
ind = length(u1)-ind+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ind == 1
    thr = estimated_gain_dist(2,ind+1);
else
    thr = estimated_gain_dist(2,ind+1);
end
if isempty(thr)
    if sum(u) == 0
        frame_length = 1;
    else
        frame_length  = ceil(sum(u));
    end
    thr = 1;
    permitted = gain >= 1;
else
    frame_length = max(ceil(sum(estimated_gain_dist(1,ind+1:end))),1);
    permitted = gain >= thr;
    
end







