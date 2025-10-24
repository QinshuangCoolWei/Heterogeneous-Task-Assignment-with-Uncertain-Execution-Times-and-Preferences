%% main
% Define parameters
% mini case
N = 4;
M = 2;  % number of members
L_vec = [1.5,1.2]; % capacity of each member
f_bar = [0.4,0.6; 0.6,0.5; 0.4,0.6; 0.6,0.7];      % f in [0, 1]
C_l = 1;
C_u = 6;
c_bar = ones(N, M) * 2.0;
c_bar(1:2,:) = 1.5;
% destined q
q = [0.35,0.3; 0.3,0.35; 0.3,0.25; 0.25,0.35];
r_bar = q .* c_bar;



alpha = 0; % approximation parameter (0: exact, >0: approximation gap)

T = 200001;

% Initialize Environment and Algorithms
% Number of trials
TT = 100;
r0 = zeros(TT, T);  % regret 
p_f = zeros(TT, T);  % penalty for f
r_sub = zeros(TT, T);  % sub-optimal regret 
% Run simulations
for s = 1:TT
    par = initialize(N,M, L_vec, C_l, C_u, T);
    [instance,feasible_assign] = generate(par, T, c_bar, r_bar,f_bar);
    [reward_opt,cumulative_rewards,cumulative_penalties_f,feedback] = simulate(par, instance, T,alpha,f_bar);
    r0(s, :) = reward_opt-cumulative_rewards;
    r_sub(s, :) = reward_opt/(1+alpha)-cumulative_rewards;
    p_f(s, :) = cumulative_penalties_f;
    q_0 = r_bar ./ c_bar;

%     figure
%     plot(r0(s, :))
%     figure
%     plot(cumulative_penalties_f)
end

% Save data
% save('multi_approx_small_50000.mat');






%% functions

function par = initialize(N,M, L_vec,C_l, C_u, T)
% initialize parameters
par.N = N;
par.M = M;
par.A_0 = find_all_action(N,M); % find all possible actions 
par.L_vec = L_vec;
par.C_u = C_u;
par.C_l = C_l;
par.T = T;
par.chosen_action = zeros(N, M);
par.release_round = zeros(N, M);
par.occupied = zeros(N, M);
par.Ts = zeros(N, M);   % # of execution
par.Tft = zeros(N, M);   % # of rounds for execution (estimate resource f)
par.q = zeros(N, M);
par.f_lcb = zeros(N,M); % LCB of f
par.Rs = zeros(N, M);
par.Ft = zeros(N, M);   % f is updated for all rounds
par.Cs = zeros(N, M);
par.Vs = zeros(N, M);
par.A_s = par.A_0;  % estimated feasible set
par.ts = 1;
par.as = [];
end


function [instance,feasible_assign] = generate(par, T, c_bar, r_bar,f_bar)
% generate random completion times, rewards, and resource occupation f
% according to avg, find corresponding best per-unit-time reward (avg)
instance.c = zeros(T, par.N, par.M);
instance.r = zeros(T, par.N, par.M);
instance.f = zeros(T, par.N, par.M);
r_range = 0.3;  % so that r in [-r_range,r_range]+r_bar
f_range = 0.15;  % so that f in [-f_range,f_range]+f_bar
for i = 1:par.N
    for j = 1:par.M
%         instance.r(:, i, j) = binornd(1, r_bar(i,j), [T, 1]);
        instance.r(:, i, j) = r_bar(i,j)-r_range + 2*r_range* rand(T, 1);
        instance.f(:, i, j) = f_bar(i,j)-f_range + 2*f_range* rand(T, 1);
        p = (c_bar(i,j) - par.C_l) / (par.C_u - par.C_l);
%         p=1;
        instance.c(:, i, j) = binornd(par.C_u - par.C_l, p, [T, 1]) + par.C_l;
    end
end
q_known = r_bar ./ c_bar;
feasible_assign = find_feasible_action(par,f_bar,0); % find true feasible actions
% find best solution given q and feasible set
best_act = opt_act(q_known,feasible_assign,par);
instance.optimal_reward = best_act.reward; % the best avg reward
instance.optimal_action = best_act.action; % the best avg reward
end


function [reward_opt,cumulative_rewards,cumulative_penalties_f,feedback] = simulate(par,instance , T, alpha,f_bar)
% find cumulative regret at each round (main function)
par = initialize(par.N,par.M, par.L_vec,par.C_l, par.C_u, T);
feedback = cell(T + par.C_u, 1);
cumulative_rewards = zeros(T, 1);
cumulative_penalties_f = zeros(T, 1);
penalty_opt = zeros(T, 1);
reward_opt = zeros(T, 1);
% regret_approx = zeros(T, 1);
cumulative_reward = 0;
cumulative_penalty_f = 0;
for t = 1:T
    action = par.chosen_action;
    resource_occupied = zeros(par.N,par.M);
    resource_occupied_avg = zeros(par.N,par.M);
    resource_opt_act = zeros(par.N,par.M);
    for i = 1:par.N
        for m = 1:par.M
            if action(i,m) == 1
                par.occupied(i,m) = 1;    % task occupied
                par.release_round(i,m) = t + instance.c(t, i,m);
                feedback{t + instance.c(t, i, m)}(i, :) = [instance.c(t, i, m), instance.r(t, i,m), t,m]; % we get the info (c,r) at t+c
            end
            if t == par.release_round(i,m)
                par.occupied(i,m) = 0;  % not performing the task
            end
            resource_opt_act(i,m) = instance.f(t,i,m)*instance.optimal_action(i,m); % resource that will be occupied by best action
            if par.occupied(i,m) == 1
                % update f_t
                par.Tft(i,m) = par.Tft(i,m) + 1;
                par.Ft(i,m) = par.Ft(i,m)+instance.f(t,i,m);
                resource_occupied(i,m) = instance.f(t,i,m); 
                resource_occupied_avg(i,m) = f_bar(i,m);                 
            end
        end
    end
%     penalty_act = sum(max(sum(resource_occupied,1) - par.L_vec,0),2);   % penalty for current occupied tasks
    penalty_act_avg = sum(max(sum(resource_occupied_avg,1) - par.L_vec,0),2);   % avg penalty for current occupied tasks
    penalty_opt_act = sum(max(sum(resource_opt_act,1) - par.L_vec,0),2);   % penalty for best assignments
    cumulative_penalty_f = cumulative_penalty_f+penalty_act_avg;
    par = update(par, t+1 ,feedback{t+1},alpha);   % update next action
    % compute suboptimal/optimal regret
    if (penalty_act_avg==0)    % count the reward only if the current assignment is feasible
        if ~isempty(feedback{t})
            cumulative_reward = cumulative_reward + sum(feedback{t}(:,2));
        end
%         reward_opt_temp = reward_opt_temp + instance.optimal_reward;
    end
    cumulative_rewards(t) = cumulative_reward;
    cumulative_penalties_f(t) = cumulative_penalty_f;
    penalty_opt(t) = penalty_opt_act; 
%     reward_opt(t) = reward_opt_temp;
    reward_opt(t) = t  * instance.optimal_reward;
    if t==T-1
        disp(par.as)
    end

end
end





function par = update(par, t, feedback,alpha)
% select action according to feedback

for i = 1:size(feedback,1)
    if feedback(i, 1) >= 1  % if this feedback is nonempty
        m = feedback(i, 4);   % record the member
        par.occupied(i,m) = 0;
        par.Ts(i,m) = par.Ts(i,m) + 1;
        par.Cs(i,m) = par.Cs(i,m) + feedback(i, 1);
        par.Vs(i,m) = par.Vs(i,m) + feedback(i, 1) ^ 2;
        par.Rs(i,m) = par.Rs(i,m) + feedback(i, 2);
    end
end

if t == par.ts+1  % phase start time

    par.q = ones(par.N, par.M) * 10000;
    for i = 1:par.N
        for m = 1:par.M
            if par.Ts(i,m) > 0    % executed times >0
                d_r = sqrt(1.5 * log(t) / max(par.Ts(i,m), 1)); %d_r
                c_hat = par.Cs(i,m) / par.Ts(i,m);
                V = par.Vs(i,m) / par.Ts(i,m) - c_hat ^ 2;  % empirical variance
                d_c = sqrt(3 * V * log(t) / max(par.Ts(i,m), 1)) + ...
                    9 * (par.C_u - par.C_l) * log(t) / max(par.Ts(i,m), 1);% d_c
                par.q(i,m) = min(1, par.Rs(i,m) / max(par.Ts(i,m), 1) + d_r) / max(par.C_l, c_hat - d_c);  % q_hat
                d_f_im = sqrt(1.5 * log(t) / max(par.Tft(i,m), 1)); % d_f_im
                par.f_lcb(i,m) = d_f_im;    % d_f_im 
            end
        end
    end
    % compute action for next phase 
    par.as = zeros(par.N,par.M);
    f_est = par.Ft ./ max(par.Tft,1);  % estimated f
    A_s_dag = find_feasible_action(par,f_est,par.f_lcb);   %  feasible set at phase s
    if ~isempty(A_s_dag)    % if there exists feasible assignments
        % choose action that returns best rewards
        par.A_s = A_s_dag;  % phase-wise feasible set
        if alpha ==0    % exact method
            assign_s = opt_act(par.q, A_s_dag, par);
            par.as = assign_s.action;
        elseif alpha == 1
            assign_s = approx_gap(par,par.q, alpha);   % compute approximated gap
            par.as = assign_s.action;
        end
    else
         % choose action that returns least infeasibility
        par.A_s = A_0;  % choose from all possible actions
        assign_s = opt_resource(par,f_est,par.f_lcb);   % \argmin_{a \in {\mc A}^0}\sum_{m} \max \{ \sum_{i} \hat{f}_{im}(t_s) a_{im} - d^f_{am}(t) ,0 \}
        par.as = assign_s.action;
    end

    if min(par.Ts,[],'all') == 0     % in initialization phase
        par.ts = par.ts + par.C_u * par.C_u + 2 * par.C_u;    % phase end time
    else
        minT = par.T;
        for i = 1:par.N
            for m = 1:par.M
                if par.as(i,m) == 1
                    minT = min(minT, par.Ts(i,m));   % find min_i T_i(y_s)
                end
            end
        end
        par.ts = par.ts + par.C_l * minT + 2 * par.C_u; % phase end time
    end
end

% action at round t
par.chosen_action = zeros(par.N, par.M);
% covered = true;

% % see if we can start new actions
% for i = 1:par.N
%     for m = 1:par.M
%         if (sum(par.occupied(i,:),2) == 1) && par.as(i,m)==0
%             covered = false;
%         end
%     end
% end

for i = 1:par.N
    if sum(par.occupied(i,:),2) == 0  &&  sum(par.as(i,:),2) == 1   %if not occupied and need to be executed
        for m = 1:par.M
            if par.as(i,m) ==1
                act_attempt = par.occupied; 
                act_attempt(i,m) = 1;   % if add this task
                 if any(all(all(par.A_s == act_attempt, 1), 2), 3)  % adding this task will be still in A_s
                    par.chosen_action(i,m) = 1; % update next action
                    par.occupied(i,m) = 1; % update occupation
                end
            end
        end
    end
end



% if covered  % we can start new actions (in as) 
% for i = 1:par.N
%     if sum(par.occupied(i,:),2) == 0  &&  sum(par.as(i,:),2) == 1   %if not occupied and need to be executed
%         src = par.as(i,:).*par.F(i,:);  % resource occupation of task i
%         if sum(member_src+src <= par.L_vec) == par.M  % adding this task will be within the src limit
%             par.chosen_action(i,:) = par.as(i,:); % update next action 
%             member_src = member_src +src; % update resource occupation
%             par.occupied(i,:) = par.as(i,:); % update occupation
%         end
%     end
% end
end
% end


%% find all "possible" action (A_0)
function A_0 = find_all_action(N,M)
% compute all possible actions with N tasks, M members

nMatrices = (M+1)^N;    % each row contains at most an 1 
A_mat = zeros(N,M,nMatrices);
count = 0;
% Generate and check all combinations using binary representation
for num = 1:2^(N*M)   % all possible combinations
    binary_num = dec2bin(num-1, N*M)-'0';   % a vector - binary representation of num
    temp_mat = reshape(binary_num,[N,M]);   % reshape into NxM matrix
    % record the matrix as feasible only if all rows has at most an 1
    if sum(sum(temp_mat,2)>1) == 0 
        count = count+1;
        A_mat(:,:,count) = temp_mat; 
    end
end
A_0 = A_mat(:,:,1:count);
end


%% find all "feasible" action (A or A_s_dagger)
function feasible_assign = find_feasible_action(par,f_est,f_lcb)
% compute the feasible set of actions based on f_est (estimated f_bar) and
% f_lcb NxM matrix used to compute LCB of f
N = par.N;
M = par.M;
A_0 = par.A_0;
L_vec = par.L_vec;
if sum(f_lcb,'all') ==0
    f_lcb = zeros(N,M);
end


N_A0 = size(A_0,3);
feas_mat = zeros(N,M,N_A0);  %pre-deposit all possible actions
count = 0;
% Generate all combinations using base-3 representation
for k = 1:N_A0   % all possible combinations
    act = A_0(:,:,k);   % temp action
    d_f_am = N*max(act.*f_lcb,[],1); % d_f_am = max_i L^2 d_f_im a_im, where we let L = N
    temp_f = sum(act.*f_est,1)-d_f_am;    % LCB of resource need from each member
    % record the matrix as feasible only if LCB of f <L_m for all column
    if sum(temp_f>L_vec)==0   
        count = count+1;
        feas_mat(:,:,count) = act; 
    end
end
feasible_assign = feas_mat(:,:,1:count);
end


%% solve for "optimal" action
function best_act = opt_act(q,feasible_assign,par)
% compute the optimal solution based on q
num_feas = size(feasible_assign,3);
% variables
reward_act = repmat(q,1,1,num_feas).*feasible_assign; 
total_reward_act = sum(sum(reward_act,1),2);    % each layer represent total avg reward for this action
[max_reward,max_reward_ind] = max(total_reward_act);    % find the best action index
best_act.action = feasible_assign(:,:,max_reward_ind); % find the best action
best_act.reward = max_reward;
end

%% when no feasible action is found, find the one with least disobey
function best_act = opt_resource(par,f_est,f_lcb)  
% \argmin_{a \in {\mc A}^0}\sum_{m} \max \{ \sum_{i} \hat{f}_{im}(t_s) a_{im} - d^f_{am}(t) ,0 \}
N = par.N;
M = par.M;
L_vec = par.L_vec;
A_0 = par.A_0;
L_vec = par.L_vec;
if sum(f_lcb,'all') ==0
    f_lcb = zeros(N,M);
end

N_A0 = size(A_0,3);
f_exceed = zeros(N_A0,1);  %pre-deposit all possible actions
% Generate all combinations using base-3 representation
for k = 1:N_A0   % all possible combinations
    act = A_0(:,:,k);   % temp action
    d_f_am = N*max(act.*f_lcb,[],1); % d_f_am = max_i L^2 d_f_im a_im, where we let L = N
    temp_f = sum(act.*f_est,1)-d_f_am;    % LCB of resource need from each member
    % record the matrix as feasible only if LCB of f <L_m for all column
    f_exceed(k) = sum(max(temp_f-L_vec, 0),2);  % all excess need of f
end

% find the best (minimun)
[min_f,min_f_ind] = min(f_exceed,[],1);
best_act.action = A_0(:,:,min_f_ind); % find the best action
best_act.penalty_f = min_f;

end
