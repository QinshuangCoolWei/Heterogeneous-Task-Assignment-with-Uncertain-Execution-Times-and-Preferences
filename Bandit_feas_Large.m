%% main

% large case
filename = 'ori_data1.xlsx';
N = 20; % # of tasks
M = 5;  % # of members
L_vec = readmatrix(filename,'Sheet','L_vec'); % capacity of each member
L = min(N,sum(L_vec));
f_bar = readmatrix(filename,'Sheet','f_bar');   % resource needed
C_l = 3;
C_u = 10;
c_bar = readmatrix(filename,'Sheet','c_bar');
r_bar = readmatrix(filename,'Sheet','r_bar');


alpha = 1; % approximation parameter (0: exact, >0: approximation gap)

T = 1000001;

% Initialize Environment and Algorithms
% Number of trials
TT = 50;
r0 = zeros(TT, T);  % regret 
p_f = zeros(TT, T);  % penalty for f
r_sub = zeros(TT, T);  % sub-optimal regret 
% Run simulations
for s = 1:TT
    par = initialize(N,M, L_vec, C_l, C_u, T);
    instance = generate(par, T, c_bar, r_bar,f_bar);
    [reward_opt,cumulative_rewards,cumulative_penalties_f,feedback] = simulate(par, instance, T,alpha,f_bar);
    r0(s, :) = reward_opt-cumulative_rewards;
    r_sub(s, :) = reward_opt/(1+alpha)-cumulative_rewards;
    p_f(s, :) = cumulative_penalties_f;
    q_0 = r_bar ./ c_bar;

end

% Save data
save('alpha_1_large_1e6.mat');





%% functions

function par = initialize(N,M, L_vec,C_l, C_u, T)
% initialize parameters
par.N = N;
par.M = M;
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
par.ts = 1;
par.as = [];
end


function [instance] = generate(par, T, c_bar, r_bar,f_bar)
% generate random completion times, rewards, and resource occupation f
% according to avg, find corresponding best per-unit-time reward (avg)
instance.c = zeros(T, par.N, par.M);
instance.r = zeros(T, par.N, par.M);
instance.f = zeros(T, par.N, par.M);
r_range = 0.15;  % so that r in [-r_range,r_range]+r_bar
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
% feasible_assign = find_feasible_action(par,f_bar,0); % find true feasible actions
% find best solution given q
best_act = opt_act(q_known,par,f_bar,0,0);
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
    par = update(par, t+1 ,feedback{t+1},alpha,T);   % update next action
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
%         disp(par.as)
    end

end
end





function par = update(par, t, feedback,alpha,T)
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
                d_c = sqrt(1.5 * V * log(t) / max(par.Ts(i,m), 1)) + ...
                    9 * (par.C_u - par.C_l) * log(t) / max(par.Ts(i,m), 1);% d_c
                par.q(i,m) = min(1, par.Rs(i,m) / max(par.Ts(i,m), 1) + d_r) / max(par.C_l, c_hat - d_c);  % q_hat
                d_f_im = sqrt(1.5 * log(t) / max(par.Tft(i,m), 1)); % d_f_im
                par.f_lcb(i,m) = d_f_im;    % d_f_im 
            end
        end
    end
    % compute action for next phase 
    par.as = zeros(par.N,par.M);
    par.f_est = par.Ft ./ max(par.Tft,1);  % estimated f
    % combine find feasible actions and find best action 
    assign_s = opt_act(par.q, par, par.f_est,par.f_lcb, alpha);   
    if ~isempty(assign_s)    % if there exists feasible assignments
        par.as = assign_s.action;
    else
         % choose action that returns least infeasibility
        assign_s = opt_resource(par, par.f_est,par.f_lcb);   % \argmin_{a \in {\mc A}^0}\sum_{m} \max \{ \sum_{i} \hat{f}_{im}(t_s) a_{im} - d^f_{am}(t) ,0 \}
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


for i = 1:par.N
    if sum(par.occupied(i,:),2) == 0  &&  sum(par.as(i,:),2) == 1   %if not occupied and need to be executed
        for m = 1:par.M
            if par.as(i,m) ==1
                act_attempt = par.occupied; 
                act_attempt(i,m) = 1;   % if add this task
                feasibility = act_in_A_s_dag(par,par.f_est,par.f_lcb, act_attempt);
                 if feasibility==1  % adding this task will be still in A_s
                    par.chosen_action(i,m) = 1; % update next action
                    par.occupied(i,m) = 1; % update occupation
                end
            end
        end
    end
end

end




%% solve for "optimal" action
function sol = opt_act(q, par,f_est,f_lcb, alpha)
% compute the optimal solution based on q, with alpha the MIP gap
N = par.N;
M = par.M;
L_vec = par.L_vec;
if sum(f_lcb,'all') ==0
    f_lcb = zeros(N,M);
end

% d_f_am = N*max(act.*f_lcb,[],1); % d_f_am = max_i L^2 d_f_im a_im, where we let L = N
% temp_f = sum(act.*f_est,1)-d_f_am;    % LCB of resource need from each member

% variables
a = binvar(N,M,'full'); % action
constraints = [ 
    sum(a.*f_est,1)-N*max(a.*f_lcb,[],1)<= L_vec,... %capability bound
    sum(a,2)<= 1     % task constraint
];

obj= -sum(q.*a,'all');  % max reward
% define MIPGap
gapdefault = 1e-4; % default value of MIP gap in gurobi
if alpha <= gapdefault
    gap = gapdefault; 
else 
    gap = alpha;
end
% ops = sdpsettings('gurobi.MIPGap',.1,'gurobi.MIPGapAbs',0.01)
options = sdpsettings('verbose',0,'solver','gurobi','gurobi.MIPGap',gap);
sol = optimize(constraints,obj,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 sol.reward = value(-obj);
 sol.action = value(a);
else
 display('Feasible solution not found');
%  sol.info
%  yalmiperror(sol.problem);
end

end

%% check action feasibility
function feasibility = act_in_A_s_dag(par,f_est,f_lcb, act)
% check whether a given action is feasible
N = par.N;
M = par.M;
L_vec = par.L_vec;
if sum(f_lcb,'all') ==0
    f_lcb = zeros(N,M);
end

feasibility = 1;
if sum(act,2) > 1     % task constraint
    feasibility =0;  
end

d_f_am = N*max(act.*f_lcb,[],1); % d_f_am = max_i L^2 d_f_im a_im, where we let L = N
temp_f = sum(act.*f_est,1)-d_f_am;    % LCB of resource need from each member
if any(temp_f >= L_vec)
    feasibility =0;
end
end


%% when no feasible action is found, find the one with least disobey
function best_act = opt_resource(par,f_est,f_lcb)  
% \argmin_{a \in {\mc A}^0}\sum_{m} \max \{ \sum_{i} \hat{f}_{im}(t_s) a_{im} - d^f_{am}(t) ,0 \}
N = par.N;
M = par.M;
L_vec = par.L_vec;
if sum(f_lcb,'all') ==0
    f_lcb = zeros(N,M);
end


% variables
a = binvar(N,M,'full'); % action
constraints = [ sum(a,2)<= 1     % task constraint
    ];
obj= sum(max(sum(act.*f_est,1)-N*max(act.*f_lcb,[],1)- L_vec, 0),'all'); % minimize exceeded resource
options = sdpsettings('verbose',0,'solver','gurobi');
sol = optimize(constraints,obj,options);

% Analyze error flags
if sol.problem == 0
 best_act.action = value(a); % find the best action
 best_act.penalty_f = value(obj);

else
 display('something wrong');
%  sol.info
%  yalmiperror(sol.problem);
end
end
