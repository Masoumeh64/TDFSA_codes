clear all
clc
global arrival_rate
global wmin

num_node = input('Please enter number of nodes:');
lambda = input('Please enter the arrival rate to each node \lambda:');
num_total_frames = input('Please enter the number of frames to simulate \n (larger number of frames result in more precise result):');
wmin = input('The minimum frame length wmin in Step 1: \n (Start from 1 and increase to get the minimum average age)')
arrival_rate = lambda*ones(num_node,1);


% The "estimated_gain_dist" is the estimated gain distribution of a
% typical node. It corresponds to the "\hat{f}^a_t" in the paper.
% estimated_gain_dist(1,i) shows the number of nodes that have age-gain
% estimated_gain_dist(2,i).
estimated_gain_dist = [0 num_node;0 1];

% The "age_AP" shows the AoI of nodes at AP at the current slot t.
% It corresponds to the "y^i_t" in the paper.
age_AP = ones(num_node,1);

% The "age_node" show the AoI of nodes at the current slot t.
% It corresponds to the "x^i_t" in the paper.
age_node = ones(num_node,1);

% Here the ages at the AP are initialized randomly.
age_AP(1:num_node) = unidrnd(num_node,num_node,1);

% The "backlogged" shows if nodes have packets at the beginning of the
% current slot.
backlogged = ones(num_node,1);

% The "col_tx" counts the number of collided slots across time.
col_tx = 0;

% The "succ_tx" counts the number of successful slots across time.
succ_tx = 0;

% The "empty_tx" counts the number of empty slots across time.
empty_tx  =0;

% The "avg_age" the average age of the nodes at AP.
avg_age = zeros(num_node,1);

% The "frame_num" counts the number of frames.
frame_num = 0;

% The t_curr shows the current time slot index.
t_curr = 100000;

% The
aa =[0 0];

% The "f_len" is a vector that stores the lenght of the last 100000 frames.
f_len =[];

% The "cnt1" counts the simulated number of slots.
% The "cnt2" counts the number of times the results are printed out.
cnt1 = 0;cnt2 = 0;

% The "w" counts the number of frames whose all slots are collided.
w = 0;

% The "n_cons_col" is a vector that stores the collided status of the
% last 20 frames. If all slots of the last frame are collided, a one
% is appended to "n_cons_col", otherwise a zero is appended.
n_cons_col = [];

% "n_reset" shows how many times we have reset the algorithm.
% The reset happens if the last 20 frames are all entirely collided.
n_reset = 0;

% The "AAR" is the sum of average age reductions across time.
% in each frame the age reduction of all nodes are divided by the
% number of nodes and the frame length and then is added to the
% previous value of AAR.
AAR = 0;

% We have a window of length 20 which slides forward frame by frame
% and in each frame, we count the number of entirely collided frames
% inside the window. If it is equal to 20, then a reset happens and a 1
% is appended to the reset vector.
reset = [];

% The "complexity" counts the number of operations.
complexity = 0;

while (1)
    
    % In this part, we check that if the last 20 frames has been entirely
    % collided, i.e., all slots of each of them are collided. If Yes,
    % then all parameters are initialized again.
    if  sum(n_cons_col) == 20
        n_reset = n_reset+1;
        
        % Here, since a reset is happend all nodes keeps only fresh
        % updates, and the "estimated_gain_dist" is reset to one for
        % all nodes
        age_node = ones(num_node,1);
        estimated_gain_dist = [0 num_node;0 1];
        backlogged = ones(num_node,1);
        avg_age = zeros(num_node,1);
        
        col_tx = 0;
        succ_tx = 0;
        empty_tx  =0;
        
        frame_num = 0;
        t_curr = 100000;
        aa =[0 0];
        f_len =[];
        cnt1 = 0;cnt2 = 0;
        w = 0;
        n_cons_col = [];
        AAR = 0;
        
        % Since we have reset the protocol, 1 is appended to "reset".
        % Also, we only keep the status of the last 10 checkings for
        % reset.
        reset = [reset 1];
        if length(reset) >10
            reset(1)= [];
        end
    else
        % Here, the 20 last frames are not all collided, i.e., we don't
        % have 20 consecutive entirely collided frames. Thus, 0 is
        % appended to the vector "reset".
        reset = [reset 0];
        if length(reset) >10
            reset(1)= [];
        end
    end
    
    
    % The "gain" keeps the age gain of the nodes.
    % It corresponds to "r^t_i" in the paper.
    gain = age_AP-age_node;
    
    % The "real_gain_dist" computes the real gain distribution of the
    % nodes through function "real_dist".
    real_gain_dist = real_dist(gain);
    
    estimated_not_tuned = estimated_gain_dist;
    
    % ****** Step 4 of the algorithm T-DFSA ***********
    % Here we truncate the estimated age-gain distribution
    % max(age_AP)-1 is the maximum possible age-gain
    % The "length(estimated_gain_dist(1,:))-1" shows the maximum
    % estimated age-gain.
    if length(estimated_gain_dist(1,:))-1 > max(age_AP)-1
        trunc = 1;
        temp = sum(estimated_gain_dist(1,1:max(age_AP)));
        estimated_gain_dist(1,1:max(age_AP)) = num_node*estimated_gain_dist(1,1:max(age_AP))/temp;
        estimated_gain_dist(:,max(age_AP)+1:end) = [];
    else
        trunc = 0;
        temp = 0;
    end
    
    %________________________________________
    % ****** Step 1 of the algorithm T-DFSA ***********
    
    % The "assign_frame_length" function outputs the frame_length
    % (corresponding to "w_t" in the paper), the threshold
    % (corresponding to "\Gamma_t" in the paper), the vector "permitted"
    % where permitted(i) shows that node i is allowed to transmit.
    if cnt1 >= 0
        [frame_length, threshold, estimated_gain_dist, permitted] = ...
            assign_frame_length(estimated_gain_dist,gain);
    end
    
    % Here, the ages at ages across the frame are summed up and added
    % to "avg_age"
    avg_age = avg_age+age_AP*frame_length+(frame_length-1)*(frame_length)/2;
    age_AP_pre = age_AP;
    
    % Here the ages at AP and nodes are updated for the end of the
    % frame. These updates are written assuming that none of the nodes
    % have been able to transmit successfully. In the next step after
    % simulating contention, the age of successful nodes will be
    % revised.
    age_AP = age_AP+frame_length;
    age_node = age_node+frame_length;
    
    
    %________________________________________________
    % Here, given "frame_length" and "threshold", the contention
    % happens.
    tx_prob = 1;
    fixed_par = [num_node, tx_prob, frame_length];
    Age = [age_AP age_node];
    
    % The "active_nodes" are those nodes who are backlogged and allowed
    % to transmit
    active_nodes = permitted.*backlogged;
    
    % The function contention returns "I" which contains the number of
    % successful, collided and empty slots in the frame, "Age" contains
    % the updated ages of AP and nodes regarding the successful
    % transmission, and "obs_gain_dist" is the observed age-gain by AP
    % (This correponds to "n^a_{t,s}" in the paper)
    [I,succ_nodes,Age,obs_gain_dist] = contention(fixed_par,Age,active_nodes);
    n_succ = I(1);n_empty = I(2); n_col = I(3);
    age_AP = Age(:,1);age_node = Age(:,2);
    backlogged(succ_nodes) = 0;
    
    %______________________________________________________
    %         % ****** Step 2 of the algorithm T-DFSA ***********
    
    % The function "optimize_num_active" receives the status of the
    % slots, i.e., "obs_slot", the threshold, and the observed age
    % gains, i.e., "obs_gain_dist", and estimates the number of
    % contending nodes with age-gain "a", in the output
    % "active_gain_dist". This corresponds to "\hat{m}^a_t" in the
    % paper.
    
    obs_slot = [n_succ, n_empty, n_col];
    active_gain_dist = optimize_num_active(fixed_par, threshold, obs_slot,obs_gain_dist);
    
    % Here, the complexity is computed
    rr = length(estimated_gain_dist)-1;
    if isempty(active_gain_dist)
        complexity = complexity+frame_length*rr;
    else
        complexity = complexity+n_succ*sum(active_gain_dist(1,:))+frame_length*rr;
    end
    
    % Here, the variables "n_cons_col" and "w" are updated.
    if n_col == frame_length
        n_cons_col = [n_cons_col 1];
        if length(n_cons_col) > 20
            n_cons_col(1) = [];
        end
        w = w+1;
    else
        n_cons_col = [n_cons_col 0];
        if length(n_cons_col) > 20
            n_cons_col(1) = [];
        end
    end
    
    %________________________________________________________________
    % Here, equations (36) and (37) of the paper are implemented, i.e.,
    % "\hat{n}^a_t" and "\hat{f}^a_t" are derived.
    
    estimated_gain_dist_plus = update_estimation(estimated_gain_dist,active_gain_dist,...
        fixed_par,obs_slot,threshold,t_curr);
    
    estimated_gain_dist_plus(1,:) = num_node*estimated_gain_dist_plus(1,:)/sum(estimated_gain_dist_plus(1,:));
    %
    %_____________________________________________________________________
    % ****** Step 3 of the algorithm T-DFSA ***********
    % Here, the function "update_estimation_vs_arrival3" derives the
    % estimated age-gain distributions considering the effect of arrival rates.
    
    estimated_gain_dist = update_estimation_vs_arrival3...
        (fixed_par,estimated_gain_dist_plus,cnt1,t_curr,obs_slot,age_AP_pre);
    
    % The function "update_backlogged_age" updates the backlogged
    % status of the nodes and their ages, considering the packet
    % generation process
    [backlogged,age_node] = update_backlogged_age(fixed_par,backlogged,age_node);
    
    if frame_length > 0
        AAR = AAR+sum(age_AP-age_AP_pre)/num_node/frame_length ;
    end
    %_____________________________________________________________________
    
    succ_tx = succ_tx+n_succ;
    col_tx = col_tx+n_col;
    empty_tx = empty_tx+n_empty;
    t_curr = t_curr+frame_length;
    frame_num = frame_num+1;
    
    throughput = succ_tx/t_curr;
    
    f_len = [f_len frame_length];
    if length(f_len) > 100000
        f_len(1) = [];
    end
    
    cnt1 = cnt1+1;
    if mod(frame_num,300) == 0
        
        % The following lines compute the distribution of frame length
        fa = unique(f_len);
        aa = repmat(f_len, length(fa),1);
        bb = repmat(fa', 1,length(f_len));
        fre = sum(aa == bb,2)';
        frame_length_dist= [fre/sum(fre);fa];
        clear aa
        clear bb
        %%%%%%%
        fprintf('Frame number: %5.3f \n', frame_num)
        fprintf('The average age (AAoI): %5.3f \n',mean(avg_age')/(t_curr-100000))
        fprintf('The N-AAoI: %5.3f \n',mean(avg_age')/(t_curr-100000)/num_node)
        fprintf('The N-AAR:  %5.3f \n',AAR/frame_num)
        fprintf('The throughput:   %5.4f \n\n',throughput)
        fprintf('----------------------------------------------------------------\n')
        
        if frame_num == num_total_frames
            break
        end
        
        
        %[complexity/frame_num complexity/(t_curr-100000)]
        
    end
    
end

