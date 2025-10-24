% Plot parameters
frame_width=3;
frame_height=2.5;
num_erb = 10;
cs = 6;
fs = 12;
R_sub_exact = '#52a447'; % green
R_sub_aprox = '#4169E1'; % blue

%% plot constrained, exact
alpha = 0.5;
T = 10000001;
num_data = 51;
r0_arr = zeros(num_data, T);  % regret 
p_f_arr = zeros(num_data, T);  % penalty for f
r_sub_arr = zeros(num_data, T);  % sub-optimal regret 
for k = 1:num_data
    filename = sprintf('large_Te7_data/Te7_%d.mat',k);
    load(filename, 'r0','r_sub','p_f');
    r0_arr(k,:) = r0;
    p_f_arr(k,:) = p_f;
    r_sub_arr(k,:) = r_sub;
    
end

% sample the whole dataset
int = 1000;

r0_arr = r0_arr(:,1:int:end);
p_f_arr = p_f_arr(:,1:int:end);
r_sub_arr = r_sub_arr(:,1:int:end);


% p_f_avg = sum(p_f_arr,1)/num_data;
% r0_avg = sum(r0_arr,1)/num_data;
% r_sub_avg = sum(r_sub_arr,1)/num_data;
%% plot infeasibility regret
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

% plot(p_f_avg,'k','LineWidth',2)


shadedErrorBar(T_arr,p_f_arr ,{@mean,@std},'lineProps','k'); 
% hold on
% shadedErrorBar(1:T,r_sub ,{@mean,@std},'lineProps','r'); 
% hold on
% shadedErrorBar(1:T,p_f ,{@mean,@std},'lineProps','g'); 
% hold off

xlabel('Round', 'FontSize', fs);
ylabel('$V_T$', 'FontSize', fs,'Interpreter','latex');
% leg = legend('$R^{inf}_T$','$R^{0.5,sub}_T$ ' ,'FontSize', fs,'Interpreter','latex','FontSize',12);
% set(leg, 'Position', [0.5, 0.1, .3, .2])
grid on;
saveas(gcf,'R_inf_avg','pdf')
saveas(gcf,'R_inf_avg','fig')


%% plot suboptimal regret

figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

% plot(r_sub_avg,'Color',R_sub_aprox,'LineWidth',2)


shadedErrorBar(T_arr,r_sub_arr ,{@mean,@std},'lineProps','b'); 

xlabel('Round', 'FontSize', fs);
ylabel('$R^{0.5}_T$', 'FontSize', fs,'Interpreter','latex');
% leg = legend('$R^{inf}_T$','$R^{0.5,sub}_T$ ' ,'FontSize', fs,'Interpreter','latex','FontSize',12);
% set(leg, 'Position', [0.5, 0.1, .3, .2])
grid on;
saveas(gcf,'large_R_sub_05_avg','pdf')
saveas(gcf,'large_R_sub_05_avg','fig')


%% plot suboptimal regret

figure;
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

plot(r0_avg)


% shadedErrorBar(1:T,r0 ,{@mean,@std},'lineProps','b'); 
% hold on
% shadedErrorBar(1:T,r_sub ,{@mean,@std},'lineProps','r'); 
% hold on
% shadedErrorBar(1:T,p_f ,{@mean,@std},'lineProps','g'); 
% hold off

xlabel('Round', 'FontSize', fs);
ylabel('$R^{0,sub}_T$', 'FontSize', fs,'Interpreter','latex');
% leg = legend('$R^{inf}_T$','$R^{0.5,sub}_T$ ' ,'FontSize', fs,'Interpreter','latex','FontSize',12);
% set(leg, 'Position', [0.5, 0.1, .3, .2])
grid on;
saveas(gcf,'R_sub_0_avg','pdf')