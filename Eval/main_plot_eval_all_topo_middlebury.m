clear all; close all;
postfix = '_topo_eval_vary_th.mat';
result_path = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/';
out_file = '/media/guoy/Research/Project_contour/MSEL_contour_extraction_full/Data/Middlebury/trainingQ_topo_eval_';
vary_p = 0.85:-0.05:0;
% eval_file1 = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB_ucm_topo_eval_vary_th.mat';
% eval_file2 = '/media/guoy/Research/Datasets/MiddEval/trainingQ_COB_MSEL_topo_eval_vary_th.mat';
% eval_file3 = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE_ucm_topo_eval_vary_th.mat';
% eval_file4 = '/media/guoy/Research/Datasets/MiddEval/trainingQ_SE_MSEL_topo_eval_vary_th.mat';
% out_file = '/media/guoy/Research/Datasets/MiddEval/trainingQ_topo_eval_all';
% 
% R_vec = [];
% P_vec = [];
% Consistency_vec = [];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(eval_file1);
% H = figure(2);
% hold on;
% plot(R, Consist, 'rs-');
% H = figure(3);
% hold on;
% plot(R, P, 'rs-');
% % plot(R_mu, P_mu, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% R_vec = [R_vec R];
% P_vec = [P_vec P];
% Consistency_vec = [Consistency_vec Consist];
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(eval_file2);
% H = figure(2);
% plot(R, Consist, 'ro--');
% H = figure(3);
% plot(R, P, 'ro--');
% % plot(R_mu, P_mu, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% R_vec = [R_vec R];
% P_vec = [P_vec P];
% Consistency_vec = [Consistency_vec Consist];
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(eval_file3);
% H = figure(2);
% plot(R, Consist, 'gs-');
% H = figure(3);
% plot(R, P, 'gs-');
% % plot(R_mu, P_mu, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
% R_vec = [R_vec R];
% P_vec = [P_vec P];
% Consistency_vec = [Consistency_vec Consist];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(eval_file4);
% H = figure(2);
% plot(R, Consist, 'go--');
% H = figure(3);
% plot(R, P, 'go--');
% % plot(R_mu, P_mu, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
% R_vec = [R_vec R];
% P_vec = [P_vec P];
% Consistency_vec = [Consistency_vec Consist];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
% H = figure(2);
% title('Consistency Evaluation')
% axis([0 1 0 1]);
% xlabel('Recall')
% ylabel('Consistency')
% grid on;
% legend('COB-ucm', 'COB-MSEL', 'SE-ucm', 'SE-MSEL', 'Location', 'SouthEast')
% set(H, 'PaperPositionMode','auto')
% out_name = [out_file '_consistency.png'];
% print(H,'-dpng','-r0', out_name);
% 
% H = figure(3);
% % plot(R_vec(1), P_vec(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% % plot(R_vec(2), P_vec(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% % plot(R_vec(3), P_vec(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
% % plot(R_mu_vec(4), P_mu_vec(4), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
% title('boundary evaluation')
% axis([0 1 0 1]);
% xlabel('Recall')
% ylabel('Precision')
% grid on;
% legend('COB-ucm', 'COB-MSEL', 'SE-ucm', 'SE-MSEL', 'Location', 'SouthEast')
% set(H, 'PaperPositionMode','auto')
% out_name = [out_file '_boundary.png'];
% print(H,'-dpng','-r0', out_name);


%%  Comparision Given SE
figure(2); hold on;
figure(3); hold on;
figure(4); hold on;
figure(5); hold on;

C_score_max = [];
C_mu = [];
C_std = [];

load([result_path 'trainingQ_SE_Kokkinos' postfix]);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'g--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'g--', 'LineWidth',2);
figure(4);
plot(R, P, 'g--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'g--', 'LineWidth',2);


load([result_path 'trainingQ_SE_Smin' postfix]);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'k--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'k--', 'LineWidth',2);
figure(4);
plot(R, P, 'k--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'k--', 'LineWidth',2);

% load([result_path 'trainingQ_SE_MinCover' postfix]);
% figure(3);
% plot(R, Consist, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);
% figure(4);
% plot(R, P, 'k--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);


load([result_path 'trainingQ_SE_MSEL' postfix]);
% load([result_path 'trainingQ_SE_MSEL_topo_eval_vary_th (another copy).mat']);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'r--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'r--', 'LineWidth',2);
figure(4);
plot(R, P, 'r--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'r--', 'LineWidth',2);

% load([result_path 'trainingQ_SE_ucm' postfix]);
load([result_path 'trainingQ_SE_ucm_topo_eval_vary_th_v2.mat']);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'c--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'c--', 'LineWidth',2);
figure(4);
plot(R, P, 'c--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'c--', 'LineWidth',2);

load([result_path 'trainingQ_SE_Kovesi' postfix]);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'b--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'b--', 'LineWidth',2);
figure(4);
plot(R, P, 'b--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'b--', 'LineWidth',2);

load([result_path 'trainingQ_SE_TCG' postfix]);
% load([result_path 'trainingQ_SE_TCG_topo_eval_vary_th (copy).mat']);
C_mu = [C_mu mean(Consist)];
C_std = [C_std std(Consist)];
figure(2);
C_scores = 2*R.*Consist./(R+Consist);
plot(vary_p, C_scores, 'm--', 'LineWidth',2);
C_score_max = [C_score_max max(C_scores)];
figure(3);
plot(R, Consist, 'm--', 'LineWidth',2);
figure(4);
plot(R, P, 'm--', 'LineWidth',2);
figure(5);
plot(C_scores, P, 'm--', 'LineWidth',2);

%% to print
H = figure(1);hold on;
hBar = bar(1:length(C_mu),C_mu);
errorbar(1:length(C_mu),C_mu, C_std,'r.')
% xt = get(gca, 'XTick');
set(gca, 'XTick', 1:length(C_mu), 'XTickLabel', {'SE-FPG', 'SE-Smin', 'SE-MSEL','SE-ucm', 'SE-Kovesi', 'SE-TCG'})
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency_hist.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);

H = figure(2);
title('Multiview Consistency')
axis([0 1 0 1]);
xlabel('Avg Boundary Probability')
ylabel('C-Measure')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency_score.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);


H = figure(3);
title('Multiview Consistency')
axis([0 1 0 1]);
xlabel('Recall')
ylabel('Consistency')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_consistency.pdf'];
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
out_name = [out_file 'SE_boundary.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);

H = figure(5);
title('Multiview Consistency')
axis([0 1 0 1]);
xlabel('C-Measure')
ylabel('Precision')
grid on;
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
% legend( 'SE-Kovesi','SE-FPG', 'SE-Smin', 'SE-MinCover', 'SE-MSEL', 'SE-ucm', 'SE-TCG', 'Location', 'SouthEast')
set(H, 'PaperPositionMode','auto')
out_name = [out_file 'SE_P_CMeasure.pdf'];
% print(H,'-dpng','-r0', out_name);
print(H,'-dpdf','-r0',out_name);
cmd = ['!pdfcrop ' out_name ' ' out_name];
eval(cmd);
