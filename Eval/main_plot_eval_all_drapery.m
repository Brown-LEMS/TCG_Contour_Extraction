clear all; close all;

eval_file1 = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB_ucm_eval.mat';
eval_file2 = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_COB_MSEL_eval.mat';
eval_file3 = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_SE_MSEL_eval.mat';
eval_file4 = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_GE_MSEL_eval.mat';
out_file = '/media/guoy/Research/Datasets/Drapery_Datasets/Cloth1_eval_all';

R_bry_mu_vec = [];
P_bry_mu_vec = [];
R_jct_mu_vec = [];
P_jct_mu_vec = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(eval_file1);
H = figure(3);
hold on;
plot(R_bry, P_bry, 'rx');
% plot(R_bry_mu, P_bry_mu, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
R_bry_mu_vec = [R_bry_mu_vec R_bry_mu];
P_bry_mu_vec = [P_bry_mu_vec P_bry_mu];

H=figure(4); 
hold on;
plot(R_jct, P_jct, 'rx');
% plot(R_jct_mu, P_jct_mu, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
R_jct_mu_vec = [R_jct_mu_vec R_jct_mu];
P_jct_mu_vec = [P_jct_mu_vec P_jct_mu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(eval_file2);
H = figure(3);
plot(R_bry, P_bry, 'gx');
% plot(R_bry_mu, P_bry_mu, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
R_bry_mu_vec = [R_bry_mu_vec R_bry_mu];
P_bry_mu_vec = [P_bry_mu_vec P_bry_mu];

H=figure(4); 
plot(R_jct, P_jct, 'gx');
% plot(R_jct_mu, P_jct_mu, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
R_jct_mu_vec = [R_jct_mu_vec R_jct_mu];
P_jct_mu_vec = [P_jct_mu_vec P_jct_mu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(eval_file3);
H = figure(3);
plot(R_bry, P_bry, 'bx');
% plot(R_bry_mu, P_bry_mu, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
R_bry_mu_vec = [R_bry_mu_vec R_bry_mu];
P_bry_mu_vec = [P_bry_mu_vec P_bry_mu];

H=figure(4); 
plot(R_jct, P_jct, 'bx');
% plot(R_jct_mu, P_jct_mu, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
R_jct_mu_vec = [R_jct_mu_vec R_jct_mu];
P_jct_mu_vec = [P_jct_mu_vec P_jct_mu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(eval_file4);
H = figure(3);
plot(R_bry, P_bry, 'mx');
% plot(R_bry_mu, P_bry_mu, 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
R_bry_mu_vec = [R_bry_mu_vec R_bry_mu];
P_bry_mu_vec = [P_bry_mu_vec P_bry_mu];

H=figure(4); 
plot(R_jct, P_jct, 'mx');
% plot(R_jct_mu, P_jct_mu, 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
R_jct_mu_vec = [R_jct_mu_vec R_jct_mu];
P_jct_mu_vec = [P_jct_mu_vec P_jct_mu];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
H = figure(3);
plot(R_bry_mu_vec(1), P_bry_mu_vec(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(R_bry_mu_vec(2), P_bry_mu_vec(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(R_bry_mu_vec(3), P_bry_mu_vec(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(R_bry_mu_vec(4), P_bry_mu_vec(4), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
title('boundary evaluation')
axis([0 1 0 1]);
xlabel('Recall')
ylabel('Precision')
grid on;
legend('COB-ucm', 'COB-MSEL', 'SE-MSEL', 'GE-MSEL', 'Location', 'NorthWest')
set(H, 'PaperPositionMode','auto')
out_name = [out_file '_bry.png'];
print(H,'-dpng','-r0', out_name);

H=figure(4); 
plot(R_jct_mu_vec(1), P_jct_mu_vec(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(R_jct_mu_vec(2), P_jct_mu_vec(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(R_jct_mu_vec(3), P_jct_mu_vec(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(R_jct_mu_vec(4), P_jct_mu_vec(4), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
title('junction evaluation')
axis([0 1 0 1]);
xlabel('Recall')
ylabel('Precision')
grid on;
legend('COB-ucm', 'COB-MSEL', 'SE-MSEL', 'GE-MSEL', 'Location', 'NorthWest')
set(H, 'PaperPositionMode','auto')
out_name = [out_file '_jct.png'];
print(H,'-dpng','-r0', out_name);