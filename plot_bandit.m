% Plot parameters
frame_width=3;
frame_height=2.5;
num_erb = 10;
cs = 6;
fs = 12;
R_sub_exact = '#52a447'; % green
R_sub_aprox = '#4169E1'; % blue

%% plot constrained, exact, p_f
% load('const_a_1_small.mat','T','r0','r_sub','p_f');
% load('const_a_0.1_small.mat','T','r0','p_f');
load('exact_small_2e6.mat','T','TT','r0','p_f')

% sample the whole dataset
int = 500;

r0_arr = r0(:,1:int:end);
p_f_arr = p_f(:,1:int:end);
% r_sub_arr = r_sub(:,1:int:end);

% p_f_avg = sum(p_f,1)/TT;
% r0_avg = sum(r0,1)/TT;
% r_sub_avg = sum(r_sub,1)/TT;
%%

N_sample = size(r0_arr,2);
T_arr = int*[0:N_sample-1]+1;

figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

shadedErrorBar(T_arr,p_f_arr ,{@mean,@std},'lineProps','k'); 

% plot(p_f_avg,'k','LineWidth',2)
xlim([-5,T]);
ylim([0,3000]);
xlabel('Round', 'FontSize', fs);
ylabel('$V_T$', 'FontSize', fs,'Interpreter','latex');

grid on;


% saveas(gcf,'small_cons_exact_pf_avg','pdf')

% plot r0
figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

shadedErrorBar(T_arr,r0_arr ,{@mean,@std},'lineProps','g'); 

% plot(r0_avg,'LineWidth',2,'Color',R_sub_exact)
xlim([-5,T]);
ylim([0,4e4]);

xlabel('Round', 'FontSize', fs);
ylabel('Regret', 'FontSize', fs,'Interpreter','latex');
grid on;


saveas(gcf,'small_cons_exact_r0_avg','pdf')





%% plot constrained, approximate
load('alpha_1_small_2e6.mat','T','TT','r0','p_f','r_sub','alpha')
% sample the whole dataset
int = 500;

r0_arr = r0(:,1:int:end);
p_f_arr = p_f(:,1:int:end);
r_sub_arr = r_sub(:,1:int:end);

% p_f_avg = sum(p_f,1)/TT;
% r0_avg = sum(r0,1)/TT;
% r_sub_avg = sum(r_sub,1)/TT;
%%

N_sample = size(r0_arr,2);
T_arr = int*[0:N_sample-1]+1;

figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);


shadedErrorBar(T_arr,p_f_arr ,{@mean,@std},'lineProps','k'); 

xlim([-5,T]);
ylim([0,3000]);
xlabel('Round', 'FontSize', fs);
ylabel('$V_T$', 'FontSize', fs,'Interpreter','latex');

grid on;


saveas(gcf,'small_cons_alpha_1_pf_avg','pdf')

%% plot r_sub
figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

% plot(r0_avg,'k','LineWidth',2, 'Color',R_sub_exact)
shadedErrorBar(T_arr,r0_arr ,{@mean,@std},'lineProps','g'); 
hold on
% plot(r_sub_avg,'k','LineWidth',2, 'Color',R_sub_aprox)
shadedErrorBar(T_arr,r_sub_arr ,{@mean,@std},'lineProps','b'); 
hold off
xlim([-5,T]);
ylim([-14e5,5e5]);
leg = legend('$R^{0}_T$','$R^{1}_T$ ' ,'FontSize', fs,'Interpreter','latex');
set(leg, 'Position', [0.2, 0.2, .3, .2])
xlabel('Round', 'FontSize', fs);
ylabel('Regret', 'FontSize', fs,'Interpreter','latex');

grid on;


saveas(gcf,'small_cons_alpha_1_r0_avg','pdf')


%% plot r1
figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);


xlim([-5,T]);
xlabel('Round', 'FontSize', fs);
ylabel('$R^{1,sub}_T$', 'FontSize', fs,'Interpreter','latex');
grid on;


saveas(gcf,'small_cons_alpha_1_r1_avg','pdf')

%%







% %% plot
% 
% alpha = 1;
% 
% load('multi_approx.mat','r0','r_sub', 'C_u','C_l','par','r_bar','c_bar','instance'); %medium Cl,Cu = 3,10
% r_1 = r0;
% r_1sub = r_sub;
% C_u1 = C_u;
% C_l1 = C_l;
% size_m = 5*20*18;
% 
% 
% T = size(r_1,2);
% T_vals = 1:T; % X-axis values
% TT = size(r0,1);
% 
% load('multi_approx2.mat','r0','r_sub', 'C_u','C_l');   %medium Cl,Cu = 1,20
% r_2 = r0;
% r_2sub = r_sub;
% C_u2 = C_u;
% C_l2 = C_l;
% 
% load('multi_approx3.mat','r0','r_sub', 'C_u','C_l');   %large
% r_3 = r0;
% r_3sub = r_sub;
% C_u3 = C_u;
% C_l3 = C_l;
% size_l = 10*20*24;
% 
% load('multi_approx_small.mat','r0','r_sub', 'C_u','C_l','par','r_bar','c_bar','instance');    %small
% r_01 = r0(:,1:size(r_3,2));
% C_u01 = C_u;
% C_l01 = C_l;
% r_01sub = r_sub;
% 
% % T = size(r_3,2);
% 
% load('multi_exact2.mat','r0', 'C'); %small
% r_0 = r0(:,1:size(r_3,2));
% C_u0 = C;
% C_l0 = 1;
% size_s = 4*2*3;
% 
% 
% q_0 = r_bar ./ c_bar;
% alg_aprox_a = approx_gap(par,q_0, alpha);
% alg_gap = sum(alg_aprox_a.*q_0,'all');
% alg_gap_T =  T_vals* (instance.optimal_reward/(1+alpha) -alg_gap);
% alg_gap_T = ones(TT,1) * alg_gap_T;
% cumulative_rewards = ones(TT,1) *T_vals* instance.optimal_reward - r_01;
% 
% shadedErrorBar(1:T,ones(TT,1)*T_vals*alg_gap-cumulative_rewards ,{@mean,@std},'lineProps','b'); 
% 
% 
% 
% % Compute mean and standard deviation for Algorithm_Proposed
% % r1_mean = mean(r_1, 1);
% % r1_std = std(r_1, 0, 1); 
% % r1_sub_mean = mean(r_1sub, 1);
% % r1_sub_std = std(r_1sub, 0, 1); 
% % 
% % r2_mean = mean(r_2, 1);
% % r2_std = std(r_2, 0, 1); 
% % r2_sub_mean = mean(r_2sub, 1);
% % r2_sub_std = std(r_2sub, 0, 1); 
% 
% %shadedErrorBar(x,2*y+20,{@mean,@std},'lineprops',{'-go','MarkerFaceColor','g'});
% % h = errorbar(T_vals, r0_mean, r0_std, 'b', 'CapSize', cs, 'DisplayName', 'Proposed');
% 
% 
% % Plot parameters
% num_erb = 10;
% cs = 6;
% fs = 12;
% % T_vals = 1:T; % X-axis values
% % transparent = 0.1;  %transparency of errorbar
% frame_width=3;
% frame_height=2.5;
% 
% 
% 
% 
% %% compare exact and approx
% 
% 
% % Plot regret for Proposed Algorithm
% figure;
% set(gcf,'Units','Inches')
% set(gcf,'Position',[4 4 frame_width frame_height])
% set(gca,'units','inches')
% set(gcf, 'PaperUnits','inches');        
% set(gcf, 'PaperSize', [frame_width frame_height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);
% % 
% % shadedErrorBar(1:T,r_01,{@mean,@std},'lineProps','g'); 
% % hold on 
% shadedErrorBar(1:T,r_01sub,{@mean,@std},'lineProps','b'); 
% hold on 
% shadedErrorBar(1:T,r_0,{@mean,@std},'lineProps','r'); 
% hold off
%  
% xlabel('Round', 'FontSize', fs);
% ylabel('Regret', 'FontSize', fs);
% leg = legend('$R^1_T$ approx','$R^0_T$ exact' ,'FontSize', fs,'Interpreter','latex','FontSize',12);
% % title('Regret Analysis for Proposed Algorithm', 'FontSize', fs);
% set(leg, 'Position', [0.3, 0.2, .3, .2])
% grid on;
% 
% % saveas(gcf,'exact_approx','pdf')
% 
% 
% 
% %% compare [C_l, C_u]
% % Plot regret for Proposed Algorithm
% figure;
% set(gcf,'Units','Inches')
% set(gcf,'Position',[4 4 frame_width frame_height])
% set(gca,'units','inches')
% set(gcf, 'PaperUnits','inches');        
% set(gcf, 'PaperSize', [frame_width frame_height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);
% 
% 
% shadedErrorBar(1:T,r_1sub,{@mean,@std},'lineProps','g'); 
% hold on
% shadedErrorBar(1:T,r_2sub,{@mean,@std},'lineProps','b'); 
% hold off
%  
% xlabel('Round', 'FontSize', fs);
% ylabel('Regret', 'FontSize', fs);
% leg = legend('$[C_l,C_u]=[3,10] $' ,'$[C_l,C_u]=[1,20]$' ,'FontSize', fs,'Interpreter','latex','FontSize',12);
% % leg.IconColumnWidth = 5;
% set(leg, 'Position', [0.42, 0.2, .2, .2])
% grid on;
% 
% % saveas(gcf,'C_l_u','pdf')
% 
% 
% %% compare size
% 
% figure;
% set(gcf,'Units','Inches')
% set(gcf,'Position',[4 4 frame_width frame_height])
% set(gca,'units','inches')
% set(gcf, 'PaperUnits','inches');        
% set(gcf, 'PaperSize', [frame_width frame_height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);
% 
% shadedErrorBar(1:T,r_01sub/8,{@mean,@std},'lineProps','g'); 
% hold on 
% shadedErrorBar(1:T,r_1sub/100,{@mean,@std},'lineProps','b'); 
% hold on 
% shadedErrorBar(1:T,r_3sub/200,{@mean,@std},'lineProps','r'); 
% hold off
%  
% xlabel('Round', 'FontSize', fs);
% ylabel('Regret', 'FontSize', fs);
% leg = legend('$N=4,M=2$' ,'$N=20,M=5$' ,'$N=20,M=10$' ,'FontSize', fs,'Interpreter','latex');
% % leg.ItemTokenSize = 0.3;
% set(leg, 'Position', [0.42, 0.2, .2, .2])
% grid on;
% 
% % saveas(gcf,'size_compare','pdf')
% 
% %% sometimes even subregret converge
% 
% figure;
% set(gcf,'Units','Inches')
% set(gcf,'Position',[4 4 frame_width frame_height])
% set(gca,'units','inches')
% set(gcf, 'PaperUnits','inches');        
% set(gcf, 'PaperSize', [frame_width frame_height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);
% 
% shadedErrorBar(1:T,r_3,{@mean,@std},'lineProps','b'); 
% hold on 
% % shadedErrorBar(1:T,r_3sub,{@mean,@std},'lineProps','r'); 
% % hold off
%  
% xlabel('Round', 'FontSize', fs);
% ylabel('Regret', 'FontSize', fs);
% leg = legend('$R^0_T$' ,'FontSize', fs,'Interpreter','latex');
% set(leg, 'Position', [0.5, 0.2, .3, .2])
% grid on;
% 
% % saveas(gcf,'R1_converge','pdf')
% 
% 



