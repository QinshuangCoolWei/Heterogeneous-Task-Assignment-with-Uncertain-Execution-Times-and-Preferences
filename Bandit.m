%% main
% Define parameters
% mini case
N = 4;
M = 2;  % number of members
L_vec = [2,1]; % capacity of each member
L = min(N,sum(L_vec));
F = ones(N,M);
F(2,1) = 1.4;
F(3,1) = 1.8; 
T = 10001;
C_l = 1;
C_u = 6;
c_bar = ones(N, M) * 5.0;
c_bar(1:2,:) = 1.5;
r_bar = ones(N, M) * 0.5;
r_bar(1,2) = 0.7;
alpha = 1; % approximation parameter (0: exact)

% larger cases
% filename = 'ori_data2.xlsx';
% N = 20;
% M = 10;  % number of members
% 
% 
% L_vec = readmatrix(filename,'Sheet','L_vec'); % capacity of each member
% L = min(N,sum(L_vec));
% F = readmatrix(filename,'Sheet','F');   % resource needed
% T = 10001;
% C_l = 1;
% C_u = 20;
% c_bar = readmatrix(filename,'Sheet','c_bar');
% r_bar = readmatrix(filename,'Sheet','r_bar');
% alpha = 1; % approximation parameter (0: exact)





% Initialize Environment and Algorithms
% Number of trials
TT = 100;
r0 = zeros(TT, T);
r_sub = zeros(TT, T);

% Run simulations
for s = 1:TT
    par = initialize(N,M, L_vec, C_l, C_u, T,F);
    instance = generate(par, T, c_bar, r_bar);
    [reward_opt,cumulative_rewards,feedback] = simulate(par, instance, T,alpha);
    r0(s, :) = reward_opt-cumulative_rewards;
    r_sub(s, :) = reward_opt/(1+alpha)-cumulative_rewards;

%     overlap=0;
%     occupation = zeros(M,T);
%     for t = 1:T
%         if ~isempty(feedback{t})
%             for i = 1:size( feedback{t},1)
%                 if feedback{t}(i,1) >=1
%                     tstart = feedback{t}(i,3);
%                     tend = feedback{t}(i,1)+tstart-1;
%                     m= feedback{t}(i,4);
%                     occupation(m,tstart:tend) = F(i,m)+occupation(m,tstart:tend);
%                 end
%             end
%            
%         end
%     end
end

% Save data
save('multi_approx_small.mat');


%% plot
load('multi_approx2.mat');

% Compute mean and standard deviation for Algorithm_Proposed
r0_mean = mean(r0, 1);
r0_std = std(r0, 0, 1); 

r1_mean = mean(r_sub, 1);
r1_std = std(r_sub, 0, 1); 


% Plot parameters
num_erb = 10;
cs = 6;
fs = 12;
T_vals = 1:T; % X-axis values
transparent = 0.1;  %transparency of errorbar

% Plot regret for Proposed Algorithm
figure;
shadedErrorBar(1:T,r_sub,{@mean,@std}); 
hold on 
shadedErrorBar(1:T,r0,{@mean,@std}); 

%shadedErrorBar(x,2*y+20,{@mean,@std},'lineprops',{'-go','MarkerFaceColor','g'});

% h = errorbar(T_vals, r0_mean, r0_std, 'b', 'CapSize', cs, 'DisplayName', 'Proposed');
xlabel('Round', 'FontSize', fs);
ylabel('Regret', 'FontSize', fs);
legend('FontSize', fs);
title('Regret Analysis for Proposed Algorithm', 'FontSize', fs);
grid on;




%% functions

function par = initialize(N,M, L_vec,C_l, C_u, T,F)
% initialize parameters
par.N = N;
par.M = M;
par.F = F;
par.L_vec = L_vec;
par.C_u = C_u;
par.C_l = C_l;
par.T = T;
par.chosen_action = zeros(N, M);
par.release_round = zeros(N, M);
par.occupied = zeros(N, M);
par.Ts = zeros(N, M);   % # of execution
par.q = zeros(N, M);
par.Rs = zeros(N, M);
par.Cs = zeros(N, M);
par.Vs = zeros(N, M);
par.ts = 1;
par.as = [];
end


function instance = generate(par, T, c_bar, r_bar)
% generate random completion times and rewards according to avg, find
% corresponding best reward (avg)
instance.c = zeros(T, par.N, par.M);
instance.r = zeros(T, par.N, par.M);
for i = 1:par.N
    for j = 1:par.M
        instance.r(:, i, j) = binornd(1, r_bar(i,j), [T, 1]);
        p = (c_bar(i,j) - par.C_l) / (par.C_u - par.C_l);
%         p=1;
        instance.c(:, i, j) = binornd(par.C_u - par.C_l, p, [T, 1]) + par.C_l;
    end
end
q_known = r_bar ./ c_bar;
% find best solution from q
sol_star = opt_gap(q_known,par);
% q_star = sort(q_star, 'descend');

instance.optimal_reward = sum(sol_star.obj); % the best avg reward
end


function [reward_opt,cumulative_rewards,feedback] = simulate(par,instance , T, alpha)
% find cumulative regret at each round (main function)
par = initialize(par.N,par.M, par.L_vec,par.C_l, par.C_u, T, par.F);
feedback = cell(T + par.C_u, 1);
cumulative_rewards = zeros(T, 1);
reward_opt = zeros(T, 1);
% regret_approx = zeros(T, 1);
cumulative_reward = 0;

for t = 1:T
    action = par.chosen_action;
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
        end
    end
    par = update(par, t+1 , feedback{t+1},alpha);   % update next action
    if ~isempty(feedback{t})
        cumulative_reward = cumulative_reward + sum(feedback{t}(:,2));
    end
    cumulative_rewards(t) = cumulative_reward;
    % compute suboptimal/optimal regret
%     regret_approx(t) = t  * instance.optimal_reward/(alpha+1) - cumulative_reward;      % sub_optimal regret
%     regret_opt(t) = t  * instance.optimal_reward - cumulative_reward;     % optimal regret
    reward_opt(t) = t  * instance.optimal_reward;
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
            end
        end
    end
    % compute action for next phase 
    par.as = zeros(par.N,par.M);
    if alpha ==0    % exact method
%         [~, sortedIndices] = sort(par.q,'descend');
%         par.as = sortedIndices(1:par.L);    % action for phase s (actions with highest q_hat)
        sol_s = opt_gap(par.q, par);
        par.as = sol_s.action;
    elseif alpha == 1
        par.as = approx_gap(par,alpha);   % compute approximated gap
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
member_src = sum(par.occupied.*par.F,1); % member resource occupied
for i = 1:par.N
    if sum(par.occupied(i,:),2) == 0  &&  sum(par.as(i,:),2) == 1   %if not occupied and need to be executed
        for m = 1:par.M
            if par.as(i,m) ==1
                src = par.as(i,m)*par.F(i,m);  % resource occupation of task i
                if member_src(m)+src <= par.L_vec(m)  % adding this task will be within the src limit
                    par.chosen_action(i,m) = 1; % update next action
                    member_src(m) = member_src(m) +src; % update resource occupation
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


%% solve for "optimal" action
function sol = opt_gap(q, par)
% compute the optimal solution based on q
F= par.F;
N = par.N;
M = par.M;
L_vec = par.L_vec;

% variables
a = binvar(N,M,'full'); % action
constraints = [ sum(a.*F,1)<= L_vec,... %capability bound
    sum(a,2)<= 1     % task constraint
];

obj= -sum(q.*a,'all');
options = sdpsettings('verbose',0,'solver','gurobi');
sol = optimize(constraints,obj,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 sol.obj = value(-obj);
 sol.action = value(a);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem);
end

end
