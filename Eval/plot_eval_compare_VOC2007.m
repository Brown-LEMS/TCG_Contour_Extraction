clear all; close all;
addpath ../util/io/

postfix = '_VOC2007_eval_vary_num_PR.mat';
% postfix_SEL = '_VOC2007_eval_vary_num_PR_len_rank.mat';
postfix_SEL = '_VOC2007_eval_vary_num_PR.mat';

%%  Comparision Given gPb
F_max = [];

figure(1); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

figure(2); hold on;
xlabel('num of curves');
ylabel('Detection Rate');
axis square;
grid on;
axis([0 500 0 1]);

load(['gPb_Kovesi' postfix]);
figure(1);
plot(Recall, Precision, 'bs-');
figure(2);
plot([vary_num 500], Recall, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['gPb_Kokkinos' postfix]);
figure(1);
plot(Recall, Precision, 'gs-');
figure(2);
plot([vary_num 500], Recall, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['gPb_KGS' postfix]);
figure(1);
plot(Recall, Precision, 'ms-');
figure(2);
plot([vary_num 500], Recall, 'ms-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['gPb_SEL' postfix_SEL]);
figure(1);
plot(Recall, Precision, 'rs-');
figure(2);
plot([vary_num 500], Recall, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(1)
legend( ['gPb-Kovesi, F-measure:' num2str(F_max(1))], ['gPb-FPG, F-measure:' num2str(F_max(2))],...
     ['gPb-KGS, F-measure:' num2str(F_max(3))], ['gPbTO-MSEL, F-measure:' num2str(F_max(4))],'Location', 'SouthEast')
figure(2)
legend( 'gPb-Kovesi','gPb-FPG', 'gPb-KGS', 'gPbTO-MSEL','Location', 'SouthEast')

%%  Comparision Given SE
F_max = [];

figure(3); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

figure(4); hold on;
xlabel('num of curves');
ylabel('Detection Rate');
axis square;
grid on;
axis([0 500 0 1]);

load(['SE_Kovesi' postfix]);
figure(3);
plot(Recall, Precision, 'bs-');
figure(4);
plot([vary_num 500], Recall, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['SE_Kokkinos' postfix]);
figure(3);
plot(Recall, Precision, 'gs-');
figure(4);
plot([vary_num 500], Recall, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['SE_SEL' postfix_SEL]);
figure(3);
plot(Recall, Precision, 'rs-');
figure(4);
plot([vary_num 500], Recall, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(3)
legend( ['SE-Kovesi, F-measure:' num2str(F_max(1))],...
     ['SE-FPG, F-measure:' num2str(F_max(2))],...
     ['SETO-MSEL, F-measure:' num2str(F_max(3))], 'Location', 'SouthEast')
figure(4)
legend( 'SE-Kovesi','SE-FPG', 'SETO-MSEL', 'Location', 'SouthEast')
%
%%  Comparision Given TO
F_max = [];
figure(5); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

figure(6); hold on;
xlabel('num of curves');
ylabel('Detection Rate');
axis square;
grid on;
axis([0 500 0 1]);


load(['TO_Kovesi' postfix]);
figure(5);
plot(Recall, Precision, 'bs-');
figure(6);
plot([vary_num 500], Recall, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_Kokkinos' postfix]);
figure(5);
plot(Recall, Precision, 'gs-');
figure(6);
plot([vary_num 500], Recall, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_SEL' postfix]);
figure(5);
plot(Recall, Precision, 'rs-');
figure(6);
plot([vary_num 500], Recall, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(5)
legend(['GTO-Kovesi, F-measure:' num2str(F_max(1))],...
    ['GTO-FPG, F-measure:' num2str(F_max(2))],...
    ['GTO-MSEL, F-measure:' num2str(F_max(3))], 'Location', 'SouthEast')
figure(6)
legend( 'GTO-Kovesi', 'GTO-FPG', 'GTO-MSEL', 'Location', 'SouthEast')


figure(1);
print_pdf('gPb_VOC2007_PR.pdf');
figure(2);
print_pdf('gPb_VOC2007_num_vs_DR.pdf');

figure(3);
print_pdf('SE_VOC2007_PR.pdf');
figure(4);
print_pdf('SE_VOC2007_num_vs_DR.pdf');

figure(5);
print_pdf('TO_VOC2007_PR.pdf');
figure(6);
print_pdf('TO_VOC2007_num_vs_DR.pdf');