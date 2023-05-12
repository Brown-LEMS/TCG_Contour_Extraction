clear all; close all;
addpath (genpath('../util/'));

prefix = 'TO_SEL_';

%% train the logistic regression for merging classifier
varying_th = (1:0.5:9)/10;
figure; hold on;
xlabel('prob th');
ylabel('error rate');
axis([0 1 0 1]);

load([prefix 'Features.mat']);
load([prefix 'Y.mat']);

Features_1 = [ones(size(Features,2),1) Features'];
Y = Y';

% Features_test = Features_1(1:length(Y)/2, :);
% Features_train = Features_1(length(Y)/2+1:end, :);
% Y_test = Y(1:length(Y)/2);
% Y_train = Y(length(Y)/2+1:end);

Features_test = Features_1;
Features_train = Features_1;
Y_test = Y;
Y_train = Y;


cnt_pos_test = sum(find(Y_test));
cnt_neg_test = sum(find(Y_test==0));

% train beta
fstd = std(Features_train);
fstd = fstd + (fstd==0);
fmean = mean(Features_train);
fmean(1) = 0;

Features_train = (Features_train - repmat(fmean,size(Features_train,1),1))./ repmat(fstd,size(Features_train,1),1);

% fit the model
fprintf(2,'Fitting model...\n');
beta = logist2(Y_train,Features_train);

beta = beta';
save_name = [prefix 'beta_of_cues_for_merging.txt'];
save(save_name,'fmean','fstd','beta','-ascii');

% test
beta = beta';
% normalize test data
Features_test = (Features_test - repmat(fmean,size(Features_test,1),1))./ repmat(fstd,size(Features_test,1),1);

p_test = 1 ./ (1 + exp(-(Features_test)*beta));  

% estimates = zeros(size(p_test));
% estimates(find(p_test>0.5)) = 1;
% 
% 
% error = sum(abs(estimates - Y_test))/length(p_test)
error = [];
for i = 1:length(varying_th)
    estimates = zeros(size(p_test));
    estimates(find(p_test> varying_th(i))) = 1;
    
    error = [error sum(abs(estimates - Y_test))/length(p_test)];
    
end

plot(varying_th, error, 'k');

min(error)

%% train using only geom continuity

fstd_1 = fstd([1 8]);
fmean_1 = fmean([1 8]);

% fit the model
fprintf(2,'Fitting model...\n');
beta_1 = logist2(Y_train,Features_train(:, [1 8]));
save_name = [prefix 'beta_of_geomcon_cue_for_merging.txt'];
beta_1 = beta_1';
save(save_name,'fmean_1','fstd_1','beta_1','-ascii');


% test
beta_1 = beta_1';

p_test_1 = 1 ./ (1 + exp(-Features_test(:, [1 8])*beta_1));  

error_1 = [];
for i = 1:length(varying_th)
    estimates = zeros(size(p_test_1));
    estimates(find(p_test_1>varying_th(i))) = 1;

    error_1 = [error_1 sum(abs(estimates - Y_test))/length(p_test_1)];
end
plot(varying_th, error_1, 'r--');

min(error_1)

%% visualize histogram
feature_names = cell(size(Features, 1), 1);
feature_names{1} = 'BgGrad';
feature_names{2} = 'SatGrad';
feature_names{3} = 'HueGrad';
feature_names{4} = 'AbsCurvature';
feature_names{5} = 'EdgeSparsity';
feature_names{6} = 'Wigg';
feature_names{7} = 'Orient';
feature_names{8} = 'Texture';

pos_features = Features(:, find(Y>0));
neg_features = Features(:, find(Y==0));

for d = 1:size(Features, 1)
    [pos_hist, pos_X] = hist(pos_features(d,:));
    [neg_hist, neg_X] = hist(neg_features(d,:));

    H = figure; hold on;
    plot(pos_X, pos_hist, '-r', 'LineWidth',8);
    plot(neg_X, neg_hist, '-b', 'LineWidth',8);
    hold off;
    AX = legend([ ' pos'], [ ' neg'], 'Location', 'NorthEast');
    LEG = findobj(AX,'type','text');
    set(LEG,'FontSize',25)
    TX = title(feature_names{d});
    TT = findobj(TX,'type','text');
    set(TT,'FontSize',25)
    print(H,'-dpng','-r0',['hist_merge/hist_' feature_names{d} '.png']);

end

P_all = 1 ./ (1 + exp(-(Features_1)*beta));
pos_P = P_all(find(Y>0));
neg_P = P_all(find(Y==0));
[pos_hist, pos_X] = hist(pos_P);
[neg_hist, neg_X] = hist(neg_P);

H = figure; hold on;
plot(pos_X, pos_hist, '-r', 'LineWidth',8);
plot(neg_X, neg_hist, '-b', 'LineWidth',8);
hold off;
AX = legend('pos', 'neg', 'Location', 'NorthEast');
LEG = findobj(AX,'type','text');
set(LEG,'FontSize',25)
TX = title('Output Prob');
TT = findobj(TX,'type','text');
set(TT,'FontSize',25)
print(H,'-dpng','-r0',['hist_merge/hist_prob' '.png']);
