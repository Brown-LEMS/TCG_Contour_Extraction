clear all; close all;
postfix = '_topo_eval_vary_num.mat';
result_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/';
out_file = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_topo_eval_';
% vary_num = 0.85:-0.05:0;


%%  Comparision Given SE
figure(2); hold on;
figure(3); hold on;
figure(4); hold on;
figure(5); hold on;

C_score_max = [];
C_mu = [];
C_std = [];

load([result_path 'trainingQ_SE_Kokkinos' postfix]);
vary_num = [vary_num 1000];
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'g--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'g--', 'LineWidth',2);
figure(4);
plot(R, P, 'g--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'g--', 'LineWidth',2);

load([result_path 'trainingQ_SE_Smin' postfix]);
vary_num = [vary_num 1000];
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'k--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'k--', 'LineWidth',2);
figure(4);
plot(R, P, 'k--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'k--', 'LineWidth',2);

% load([result_path 'trainingQ_SE_MinCover' postfix]);
% figure(3);
% plot(R, Consist, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);
% figure(4);
% plot(R, P, 'k--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);


load([result_path 'trainingQ_SE_MSEL' postfix]);
vary_num = [vary_num 1000];
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'r--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'r--', 'LineWidth',2);
figure(4);
plot(R, P, 'r--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'r--', 'LineWidth',2);

% load([result_path 'trainingQ_SE_ucm' postfix]);
load([result_path 'trainingQ_SE_ucm_topo_eval_vary_th_v2.mat']);
% vary_num = vary_num;
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'c--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'c--', 'LineWidth',2);
figure(4);
plot(R, P, 'c--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'c--', 'LineWidth',2);

load([result_path 'trainingQ_SE_Kovesi' postfix]);
vary_num = [vary_num 1000];
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'b--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'b--', 'LineWidth',2);
figure(4);
plot(R, P, 'b--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'b--', 'LineWidth',2);

load([result_path 'trainingQ_SE_TCG' postfix]);
vary_num = [vary_num 1000];
Consist_all = cntConsist_bry(:)./sumConsist_bry(:);
C_mu = [C_mu mean(Consist_all)];
C_std = [C_std std(Consist_all)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_num, C_scores, 'm--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(vary_num, Consist, 'm--', 'LineWidth',2);
figure(4);
plot(R, P, 'm--', 'LineWidth',2);
figure(5);
plot(vary_num, R, 'm--', 'LineWidth',2);

%% to print
H = figure(1);hold on;
hBar = bar(1:length(C_mu),C_mu);
errorbar(1:length(C_mu),C_mu, C_std,'r.')
% xt = get(gca, 'XTick');
set(gca, 'XTick', 1:length(C_mu), 'XTickLabel', {'SE-FPG', 'SE-Smin','SE-MSEL','SE-ucm', 'SE-Kovesi', 'SE-TCG'})
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency_hist_vary_num.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);

H = figure(2);
title('Multiview Consistency')
axis([0 1000 0 1]);
xlabel('N contours')
ylabel('C-Measure')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency_score_vary_num.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);


H = figure(3);
title('Multiview Consistency')
axis([0 1000 0 1]);
xlabel('N contours')
ylabel('Consistency')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency_vary_num.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);


H = figure(4);
title('Occlusion Boundary Identification')
axis([0 1 0 1]);
xlabel('Recall')
ylabel('Precision')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_boundary_vary_num.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);

H = figure(5);
title('Recover of occlusion boundary')
axis([0 1000 0 1]);
xlabel('N contours')
ylabel('Recall')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_Recall_N_vary_num.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);