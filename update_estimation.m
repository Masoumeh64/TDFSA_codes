function estimated_gain_dist = update_estimation(estimated_gain_dist,active_gain_dist,...
    fixed_par,obs_slot, thr,t_curr,age_AP)
global arrival_rate
lam = arrival_rate(1);
n_succ = obs_slot(1);
n_col = obs_slot(3);
num_node = fixed_par(1);
frame_length = fixed_par(3);
% estimated_gain_dist
if isempty(active_gain_dist)
    estimated_gain_dist(1,thr+1:end) = 0;
else
    
    u = length(estimated_gain_dist(2,thr+1:end));
    estimated_gain_dist(1,thr+1:end) = zeros(1,u);
    ind = active_gain_dist(2,:)+1;
    
    
%     if n_col ~= frame_length
        estimated_gain_dist(1,ind) = active_gain_dist(1,:);
        estimated_gain_dist(2,ind(1):ind(end)) = (ind(1):ind(end))-1;
        estimated_gain_dist(2,1:ind(1)-1) = 0:ind(1)-2;
%     else
%         age_max = max(age_AP)-1;
%         n = 1-(1-lam)^(t_curr+2-thr);
%         if n == 0
%             1;
%         end
%         estimated_gain_dist(1,thr+1:t_curr+2) = ...
%             active_gain_dist(1,:)*lam*(1-lam).^(0:t_curr-thr+1)/n;
%         estimated_gain_dist(2,thr+1:t_curr+2) = thr:t_curr+1;
%     end
    
end
estimated_gain_dist(1,1) = estimated_gain_dist(1,1)+n_succ;

ind1 = find(estimated_gain_dist(2,2:end)==0);
if ~isempty(ind1)
    1;
end


%_________________________________________________________________________
% no_new_gains = size(active_gain_dist,2);
%
% for i = 1:no_new_gains
%     if isnan(max_gain)
%         estimated_gain_dist = active_gain_dist;
%     else
%         if active_gain_dist(2,i) > max_gain
%             estimated_gain_dist = [estimated_gain_dist active_gain_dist(:,i)];
%         elseif gain_dist(2,i) < min_gain
%             estimated_gain_dist = [active_gain_dist(:,i) estimated_gain_dist];
%         else
%             ind = find(estimated_gain_dist(2,:) == active_gain_dist(2,i));
%             estimated_gain_dist(1,ind) = estimated_gain_dist(1,ind)+1;
%         end
%     end
% end
%
% if size(estimated_gain_dist,2)~= 0 && estimated_gain_dist(1,end) == 0
%     estimated_gain_dist(:,end) = [];
% end