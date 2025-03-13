function as = approx_gap(par, alpha)   
 % compute approximated gap using methods from "An efficient approximation
 % for the Generalized Assignment Problem" (alg 2)
F= par.F;
N = par.N;
M = par.M;
L_vec = par.L_vec;
q = par.q;  % profit p

P = zeros(N,M);
S_j_bar = zeros(N,M);

T_bin = -ones(N,1); % initialization
% 2.a
for m = 1:M
    for i = 1:N
        if T_bin(i) == -1
            P(i,m) = q(i,m);
        else 
            P(i,m) = q(i,m)-q(i,T_bin(i)); % if T_bin(i) =k
        end
    end
    % 2.b run knapsack alg on bin m
    if alpha ==1
        f = -P(:,m); % coefficient for obj fcn (min)
        intcon = 1:N; % integers count
        A = [F(:,m)]';
        b = L_vec(m);
        lb = zeros(N,1);    % specify binary bounds
        ub = ones(N,1);
        Aeq=[];
        beq=[];
        x0=[];
        options = optimoptions('intlinprog','Display','off');
        x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x0,options);  % find knapsack sol
        S_j_bar(:,m) = x;
%         S_j_bar(:,m) = knapsack_exact(P(:,m),L_vec(m),F(:,m)); % S_j_bar(i,m) = 1 if assigned, 0 otherwise
    else
    end

    
    % 2.c
    T_bin(S_j_bar(:,m)==1)=m;    % T_bin(i)=m if  S_j_bar(i,m)=1
end

% 3. final assignment i is mapped to T_bin(i) if ～= -1; 
as =  zeros(N,M);
for m = 1:M
    as(:,m)= double(T_bin == m);
end
end

%%
% function assignment = knapsack_exact(par)
% F= par.F;
% N = par.N;
% M = par.M;
% L_vec = par.L_vec;
% q = par.q;  % profit p
% 
% end