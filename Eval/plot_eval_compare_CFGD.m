clear all; close all;
addpath ../util/io/
postfix = '_CFGD_eval_vary_num_PR.mat';

% %% comparision of using TO and not TO
% figure(1); hold on;
% xlabel('Recall');
% ylabel('Precision');
% axis square;
% grid on;
% axis([0 1 0 1]);
% 
% figure(2); hold on;
% xlabel('num of curves');
% ylabel('Recall');
% axis square;
% grid on;
% axis([0 500 0 1]);
% 
% load(['gPb_SEL' postfix]);
% figure(1);
% plot(Recall, Precision, 'rs-');
% figure(2);
% plot([vary_num 500], Recall, 'rs-');
% 
% load(['gPb_NO_TO_SEL' postfix]);
% figure(1);
% plot(Recall, Precision, 'ro-');
% figure(2);
% plot([vary_num 500], Recall, 'ro-');
% 
% load(['gPb_Kokkinos' postfix]);
% figure(1);
% plot(Recall, Precision, 'gs-');
% figure(2);
% plot([vary_num 500], Recall, 'gs-');
% 
% load(['gPb_Kovesi' postfix]);
% figure(1);
% plot(Recall, Precision, 'bs-');
% figure(2);
% plot([vary_num 500], Recall, 'bs-');
% 
% load(['gPb_KGS' postfix]);
% figure(1);
% plot(Recall, Precision, 'ms-');
% figure(2);
% plot([vary_num 500], Recall, 'ms-');
% 
% figure(1)
% legend('gPb-TO-SEL', 'gPb-SEL', 'gPb-FPG', 'gPb-Kovesi', 'gPb-KGS', 'Location', 'SouthEast')
% figure(2)
% legend('gPb-TO-SEL', 'gPb-SEL', 'gPb-FPG', 'gPb-Kovesi', 'gPb-KGS', 'Location', 'SouthEast')
% 
% keyboard;
% close all;

% %%  Comparision Given gPb
% F_max = [];
% 
% figure(1); hold on;
% xlabel('Recall');
% ylabel('Precision');
% axis square;
% grid on;
% axis([0 1 0 1]);
% 
% figure(2); hold on;
% xlabel('num of curves');
% ylabel('Recall');
% axis square;
% grid on;
% axis([0 500 0 1]);
% 
% load(['gPb_Kovesi' postfix]);
% figure(1);
% plot(Recall, Precision, 'bs-');
% figure(2);
% plot([vary_num 500], Recall, 'bs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% load(['gPb_Kokkinos' postfix]);
% figure(1);
% plot(Recall, Precision, 'gs-');
% figure(2);
% plot([vary_num 500], Recall, 'gs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% load(['gPb_KGS' postfix]);
% figure(1);
% plot(Recall, Precision, 'ms-');
% figure(2);
% plot([vary_num 500], Recall, 'ms-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% load(['gPb_SEL' postfix]);
% figure(1);
% plot(Recall, Precision, 'rs-');
% figure(2);
% plot([vary_num 500], Recall, 'rs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% figure(1)
% legend(['gPb-Kovesi, F-measure:' num2str(F_max(1))],...
%     ['gPb-FPG, F-measure:' num2str(F_max(2))],...
%     ['gPb-KGS, F-measure:' num2str(F_max(3))], ...
%     ['gPbTO-MSEL, F-measure:' num2str(F_max(4))], 'Location', 'SouthEast')
% figure(2)
% legend('gPb-Kovesi', 'gPb-FPG',  'gPb-KGS', 'gPbTO-MSEL','Location', 'SouthEast')


%%  Comparision Given SE
F_max = [];

figure(1); hold on;
xlabel('N contours');
ylabel('F-measure');
axis square;
grid on;
axis([0 500 0 1]);

figure(2); hold on;
xlabel('Recall');
ylabel('F-measure');
axis square;
grid on;
axis([0 1 0 1]);

figure(3); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

figure(4); hold on;
xlabel('N contours');
ylabel('Recall');
axis square;
grid on;
axis([0 500 0 1]);

load(['SE_Kokkinos' postfix]);
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'g--', 'LineWidth',2);
figure(4);
plot([vary_num 1000], Recall, 'g--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'g--', 'LineWidth',2);
figure(1);
plot([vary_num 1000], F, 'g--', 'LineWidth',2);

load(['SE_Smin' postfix]);
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'k--', 'LineWidth',2);
figure(4);
plot([vary_num 1000], Recall, 'k--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'k--', 'LineWidth',2);
figure(1);
plot([vary_num 1000], F, 'k--', 'LineWidth',2);

% load(['SE_MinCover' postfix]);
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% figure(3);
% plot(Recall, Precision, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);
% figure(4);
% plot([vary_num 1000], Recall, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);
% figure(2);
% plot(Recall, F, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);
% figure(1);
% plot([vary_num 1000], F, '--', 'color', [0.8500, 0.3250, 0.0980], 'LineWidth',2);

% load(['SE_SEL' postfix]);
load(['SE_MSEL' postfix]);
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'r--', 'LineWidth',2);
figure(4);
plot([vary_num 1000], Recall, 'r--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'r--', 'LineWidth',2);
figure(1);
plot([vary_num 1000], F, 'r--', 'LineWidth',2);

% load(['SE_ucm' postfix]);
load('SE_ucm_CFGD_eval_vary_prob_PR_v2.mat');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'c--', 'LineWidth',2);
figure(4);
plot(vary_num, Recall, 'c--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'c--', 'LineWidth',2);
figure(1);
plot(vary_num, F, 'c--', 'LineWidth',2);

load(['SE_Kovesi' postfix]);
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'b--', 'LineWidth',2);
figure(4);
plot([vary_num 1000], Recall, 'b--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'b--', 'LineWidth',2);
figure(1);
plot([vary_num 1000], F, 'b--', 'LineWidth',2);

load(['SE_TCG' postfix]);
% load('SE_TCG_CFGD_eval_vary_num_PR (copy).mat')
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];
figure(3);
plot(Recall, Precision, 'm--', 'LineWidth',2);
figure(4);
plot([vary_num 1000], Recall, 'm--', 'LineWidth',2);
figure(2);
plot(Recall, F, 'm--', 'LineWidth',2);
figure(1);
plot([vary_num 1000], F, 'm--', 'LineWidth',2);

figure(1)
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
figure(2)
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
figure(3)
legend(['Kokkinos, F-measure:' num2str(F_max(1))],...
    ['Leordeanu et al., F-measure:' num2str(F_max(2))],...
    ['Guo et al., F-measure:' num2str(F_max(3))],...
    ['UCM, F-measure:' num2str(F_max(4))],...
    ['Kovesi, F-measure:' num2str(F_max(5))],...
    ['TCG (proposed), F-measure:' num2str(F_max(6))], 'Location', 'SouthWest')
figure(4)
legend( 'Kokkinos', 'Leordeanu et al.', 'Guo et al.', 'UCM', 'Kovesi', 'TCG (proposed)', 'Location', 'SouthEast')
%

% %%  Comparision Given TO
% F_max = [];
% figure(5); hold on;
% xlabel('Recall');
% ylabel('Precision');
% axis square;
% grid on;
% axis([0 1 0 1]);
% 
% figure(6); hold on;
% xlabel('num of curves');
% ylabel('Recall');
% axis square;
% grid on;
% axis([0 1000 0 1]);
% 
% load(['TO_Kovesi' postfix]);
% figure(5);
% plot(Recall, Precision, 'bs-');
% figure(6);
% plot([vary_num 1000], Recall, 'bs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% load(['TO_Kokkinos' postfix]);
% figure(5);
% plot(Recall, Precision, 'gs-');
% figure(6);
% plot([vary_num 1000], Recall, 'gs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% load(['TO_SEL' postfix]);
% figure(5);
% plot(Recall, Precision, 'rs-');
% figure(6);
% plot([vary_num 1000], Recall, 'rs-');
% F = Recall.*Precision./(Recall + Precision)*2;
% F_max = [F_max max(F)];
% 
% figure(5)
% legend(['GTO-Kovesi, F-measure:' num2str(F_max(1))],...
%     ['GTO-FPG, F-measure:' num2str(F_max(2))],...
%     ['GTO-MSEL, F-measure:' num2str(F_max(3))], 'Location', 'SouthEast')
% figure(6)
% legend( 'GTO-Kovesi', 'GTO-FPG', 'GTO-MSEL', 'Location', 'SouthEast')
% %

% H=figure(1);
% out_file = 'gPb_CFGD_PR.pdf';
% print(H,'-dpdf','-r0',out_file);
% cmd = ['!pdfcrop ' out_file ' ' out_file];
% eval(cmd);
% % print_pdf('gPb_CFGD_PR.pdf');
% H=figure(2);
% out_file = 'gPb_CFGD_num_vs_DR.pdf';
% print(H,'-dpdf','-r0',out_file);
% cmd = ['!pdfcrop ' out_file ' ' out_file];
% eval(cmd);
% % print_pdf('gPb_CFGD_num_vs_DR.pdf');

H=figure(1);
out_file = 'SE_CFGD_F_num.pdf';
print(H,'-dpdf','-r0',out_file);
cmd = ['!pdfcrop ' out_file ' ' out_file];
eval(cmd);
H=figure(2);
out_file = 'SE_CFGD_F_R.pdf';
print(H,'-dpdf','-r0',out_file);
cmd = ['!pdfcrop ' out_file ' ' out_file];
eval(cmd);
H=figure(3);
out_file = 'SE_CFGD_PR.pdf';
print(H,'-dpdf','-r0',out_file);
cmd = ['!pdfcrop ' out_file ' ' out_file];
eval(cmd);
% print_pdf('SE_CFGD_PR.pdf');
H=figure(4);
out_file = 'SE_CFGD_num_vs_DR.pdf';
print(H,'-dpdf','-r0',out_file);
cmd = ['!pdfcrop ' out_file ' ' out_file];
eval(cmd);
% print_pdf('SE_CFGD_num_vs_DR.pdf');

% H=figure(5);
% out_file = 'TO_CFGD_PR.pdf';
% print(H,'-dpdf','-r0',out_file);
% cmd = ['!pdfcrop ' out_file ' ' out_file];
% eval(cmd);
% % print_pdf('TO_CFGD_PR.pdf');
% H=figure(6);
% out_file = 'TO_CFGD_num_vs_DR.pdf';
% print(H,'-dpdf','-r0',out_file);
% cmd = ['!pdfcrop ' out_file ' ' out_file];
% eval(cmd);
% % print_pdf('TO_CFGD_num_vs_DR.pdf');
