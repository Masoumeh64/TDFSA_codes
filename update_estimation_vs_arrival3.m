function estimated_gain_dist = update_estimation_vs_arrival3...
    (fixed_par,estimated_gain_dist,cnt1,t_curr,obs_slot,age_AP)
global arrival_rate
global wmin

lam = arrival_rate(1);
num_node = fixed_par(1);
tx_prob = fixed_par(2);
frame_length = fixed_par(3);
max_age_AP = max(age_AP);

L = estimated_gain_dist(1,:);
f_len = frame_length;
v = [];
r_pre = 0;
l = 1;

lam_vec = (1-lam).^(0:1:max_age_AP+f_len);
hmax_vec = min(t_curr+1,max_age_AP-(0:max_age_AP-1));
SS = sum(max_age_AP-(1:max_age_AP) > t_curr+1);
Hb_vec = 1-(1-lam).^hmax_vec;
D = (lam^2)*L./Hb_vec(1:length(L))/((1-lam)^2);
L_org = L;
j=0;
r = 0;

if SS > 0
    while (1)
        if length(L) < j+1
            L = [L zeros(1,j+1-length(L))];
        end
        m = repmat(1:f_len,min(j-1,max_age_AP-1)+1,1);
        b = repmat((0:min(j-1,max_age_AP-1))',1,f_len);
        h_max = min(t_curr+1,max_age_AP-b);
        m_low  = max(1,b-j+f_len+1);
        m_up = min(f_len,h_max+f_len+b-j);
        
        m1 = m <= m_up;
        m2 = m >= m_low;
        
        h   = j-b+m-f_len;
        gg  = 1-(1-lam).^h_max;
        p_h = (lam*((1-lam).^(h-1))./gg).*m1.*m2;
        p_h(isnan(p_h)) = 0;
        p_m = lam.*((1-lam).^(m-1)).*m1.*m2;
        H = p_m.*p_h;
        v(j+1) = ((1-lam)^f_len)*L(j+1)+L(1:min(j-1,max_age_AP-1)+1)*sum(H,2);
        
        [length(v)-1 max(age_AP)+f_len-1];
        if length(v)-1 > max(age_AP)+f_len-1
            1;
        end
        r_pre = r;
        r = sum(v);
        if j >= t_curr+f_len+2
            1;
        end
        j = j+1;
        if num_node-r < 1
            v(end) = num_node-sum(v(1:end-1));
            break
        end
        
    end
    
end
if SS == 0
    v1 = v;
    v = 0;
end

if SS == 0
    
    L = L_org;
    j = 1;
    com1 = 0;
    com2 = 0;
    while (1)
        if length(L) < j
            L = [L zeros(1,j-length(L))];
            D = [D zeros(1,j-length(D))];
        end
        if length(L) ~= length(D)
            1;
        end
        if j <= f_len+1
            
            y1 = 1:j-2;
            y2 = y1-(j-1)+f_len+2;
            
            if j > 1
                temp = D(j-1)*(1-lam)^(f_len+1)+D(y1)*lam_vec(y2)'+(1-lam)*com1;
            else
                temp = 0;
            end
            v(j) = (1-lam)^f_len*L(j)+temp;
            com1 = temp;
        elseif j > f_len+1 && j <= max_age_AP+1
            
            if j == f_len+2
                y1 = 1:f_len;
                y2 = 2*y1+(j-1)-f_len+1;
                y3 = 2*y1+(j-1)-1-f_len+1;
                
                temp1 = D(1)*sum(lam_vec(y2))+D(2)*sum(lam_vec(y3));
                temp2 = 0;
                for b = 2:f_len
                    y4 = 2*(b:f_len)+(j-1)-b-f_len+1;
                    temp2 = temp2+D(b+1)*sum(lam_vec(y4));
                end
                v(j) = (1-lam)^f_len*L(j)+temp1+temp2;
                com1 = temp1;
                com2 = temp2;
            else
                
                y1 = 1:f_len;
                y2 = j-f_len:j-2;
                y3 = y2-(j-1)+f_len+2;
                
                temp1 = D(j-f_len)*sum(lam_vec(2*y1+1))+(1-lam)*com1;
                temp2 = D(j-1)*(1-lam)^(f_len+1)-D(j-f_len)*sum(lam_vec(2*(y1)+1))+...
                    D(y2)*lam_vec(y3)'+(1-lam)*com2;
                
                v(j) = (1-lam)^f_len*L(j)+temp1+temp2;
                com1 = temp1;
                com2 = temp2;
            end
            
        elseif j > max_age_AP+1 && j <= max_age_AP+f_len
            if j == max_age_AP+2
                y1 = 1:f_len-1;
                temp1 = 0;
                for b = 0: max_age_AP+1-f_len
                    y2 = 2*y1+max_age_AP+1-b-f_len+1;
                    temp1 = temp1+D(b+1)*sum(lam_vec(y2));
                end
                temp2 = 0;
                for b = max_age_AP-f_len+2:max_age_AP-1
                    y1 = b-max_age_AP+f_len:f_len-1;
                    y2 = 2*y1+max_age_AP+1-b-f_len;
                    temp2 = temp2+D(b+1)*sum(lam_vec(y2+1));
                end
                v(j) = (1-lam)^f_len*L(j)+temp1+temp2;
                com1 = temp1;
                com2 = temp2;
                
            else
                y1 = 1:max_age_AP-(j-1)+f_len;
                y2 = 1:j-1-f_len;
                y3 = 2*max_age_AP-(j-1)+2+f_len-(y2-1);
                
                temp1 = D(j-f_len)*sum(lam_vec(2*y1+1))-D(y2)*lam_vec(y3+1)'+(1-lam)*com1;
                
                y1 = 1:max_age_AP-(j-1)+f_len;
                y2 = j-f_len:max_age_AP;
                y3 = y2-j+f_len+2;
                y4 = 2*max_age_AP-(j-1)+2-(y2-1)+f_len+1;
                
                temp2 = D(y2)*lam_vec(y3+1)'-D(y2)*lam_vec(y4)'...
                    -D(j-f_len)*sum(lam_vec(2*y1+1))+(1-lam)*com2;
                
                v(j) = (1-lam)^f_len*L(j)+temp1+temp2;
                
                com1 = temp1;
                com2 = temp2;
            end
        end
        r_pre = r;
        r = sum(v);
        j = j+1;
        
        if num_node-r < wmin-0.01
            v(end) = num_node-sum(v(1:end-1));
            break
        end
    end
end
estimated_gain_dist = [v;0:length(v)-1];


