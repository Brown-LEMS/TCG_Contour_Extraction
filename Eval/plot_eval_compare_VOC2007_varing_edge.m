clear all; close all;
addpath ../util/io/
postfix = '_VOC2007_eval_vary_num_PR.mat';

%%  Comparision Given MSEL varing edges
F_max = [];

figure(1); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);


load(['gPb_SEL' postfix]);
figure(1);
plot(Recall, Precision, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['SE_SEL' postfix]);
figure(1);
plot(Recall, Precision, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_SEL' postfix]);
figure(1);
plot(Recall, Precision, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(1)
legend(['gPbTO-MSEL, F-measure:' num2str(F_max(1))],...
    ['SETO-MSEL, F-measure:' num2str(F_max(2))],...
    ['GTO-MSEL, F-measure:' num2str(F_max(3))], ...
    'Location', 'SouthEast')

%%  Comparision Given Kovesi varing edges
F_max = [];

figure(2); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);


load(['gPb_Kovesi' postfix]);
figure(2);
plot(Recall, Precision, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['SE_Kovesi' postfix]);
figure(2);
plot(Recall, Precision, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_Kovesi' postfix]);
figure(2);
plot(Recall, Precision, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(2)
legend(['gPbTO-Kovesi, F-measure:' num2str(F_max(1))],...
    ['SETO-Kovesi, F-measure:' num2str(F_max(2))],...
    ['GTO-Kovesi, F-measure:' num2str(F_max(3))], ...
    'Location', 'SouthEast')

%%  Comparision Given Kovesi varing edges
F_max = [];

figure(3); hold on;
xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);


load(['gPb_Kokkinos' postfix]);
figure(3);
plot(Recall, Precision, 'bs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['SE_Kokkinos' postfix]);
figure(3);
plot(Recall, Precision, 'gs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

load(['TO_Kokkinos' postfix]);
figure(3);
plot(Recall, Precision, 'rs-');
F = Recall.*Precision./(Recall + Precision)*2;
F_max = [F_max max(F)];

figure(3)
legend(['gPbTO-FPG, F-measure:' num2str(F_max(1))],...
    ['SETO-FPG, F-measure:' num2str(F_max(2))],...
    ['GTO-FPG, F-measure:' num2str(F_max(3))], ...
    'Location', 'SouthEast')
%


figure(1);
print_pdf('MSEL_VOC2007_PR.pdf');
figure(2);
print_pdf('Kovesi_VOC2007_PR.pdf');
figure(3);
print_pdf('FPG_VOC2007_PR.pdf');
