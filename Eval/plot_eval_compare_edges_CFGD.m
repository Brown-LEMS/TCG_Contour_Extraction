clear all; close all;
addpath ../util/io/
postfix = '_CFGD_eval_vary_num_PR.mat';


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
ylabel('Recall');
axis square;
grid on;
axis([0 1000 0 1]);

load(['gen_no_interp_SEL' postfix]);
figure(5);
plot(Recall, Precision, 'bs-');
figure(6);
plot([vary_num 1000], Recall, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['gen_interp_SEL' postfix]);
figure(5);
plot(Recall, Precision, 'gs-');
figure(6);
plot([vary_num 1000], Recall, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_SEL' postfix]);
figure(5);
plot(Recall, Precision, 'rs-');
figure(6);
plot([vary_num 1000], Recall, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(5)
legend(['gradient(no interpolated)-MSEL, F-measure:' num2str(F_max(1))],...
    ['gradient(interpolated)-MSEL, F-measure:' num2str(F_max(2))],...
    ['GTO-MSEL, F-measure:' num2str(F_max(3))], 'Location', 'SouthEast')
figure(6)
legend( 'gradient(no interpolated)-MSEL', 'gradient(interpolated)-MSEL', 'GTO-MSEL', 'Location', 'SouthEast')
%

figure(5);
print_pdf('MSEL_CFGD_vary_edge_detector_PR.pdf');
figure(6);
print_pdf('MSEL_CFGD_vary_edge_detector_num_vs_DR.pdf');
