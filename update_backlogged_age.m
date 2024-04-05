function [backlogged,age_node] = update_backlogged_age(fixed_par,backlogged,age_node)
global arrival_rate

num_node = fixed_par(1);
tx_prob = fixed_par(2);
frame_length = fixed_par(3);

pack_gen = rand(num_node,frame_length) < repmat(arrival_rate,1,frame_length);

h5 = sum(pack_gen,2) > 0;
h6 = h5.*(1:num_node)';
nodes_new_pack = h6(h6>0);

h7 = pack_gen(nodes_new_pack,:).*...
    (repmat((frame_length:-1:1),length(nodes_new_pack),1));
h8 = min((h7./h7).*h7,[],2,'omitnan');

age_node(nodes_new_pack) = h8;
backlogged(nodes_new_pack) = 1;


