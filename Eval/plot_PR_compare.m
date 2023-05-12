clear all; %close all;


figure; hold on;

load('TO_SEL_CFGD_loose_eval_PR_mat.mat');
Recall = PR_mat(1,:);
Precision = PR_mat(2,:);
plot(Recall, Precision, 'rs-');

load('gPb_SEL_CFGD_loose_eval_PR_mat.mat');
Recall = PR_mat(1,:);
Precision = PR_mat(2,:);
plot(Recall, Precision, 'gs-');

load('SE_SEL_CFGD_loose_eval_PR_mat.mat');
Recall = PR_mat(1,:);
Precision = PR_mat(2,:);
plot(Recall, Precision, 'bs-');

hold off;

xlabel('Recall');
ylabel('Precision');
axis square;
grid on;
axis([0 1 0 1]);

legend('TO-SEL', 'gPb-SEL', 'SE-SEL');