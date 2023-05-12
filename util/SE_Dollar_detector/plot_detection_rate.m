clear
load Cnt-val-SE_TOrient_YL_Detection_accuracy_a_0.65_b_0.75_thr_0.7.mat
xDAEB=xs;
RDAEB = R;
load IoU-val-SE_TOrient_YL_Detection_rate_a_0.65_b_0.75.mat
RDREB = R;
xDREB=xs;
load Cnt-val-SE_TOrient_YL_Detection_accuracy_a_0.65_b_0.75_thr_0.7.mat
xDAYL=xs;
RDAYL = R;
load IoU-val-SE_TOrient_YL_Detection_rate_a_0.65_b_0.75.mat
RDRYL = R;
xDRYL=xs;
load Cnt-val-SE_TOrient_YL_Detection_accuracy_a_0.65_b_0.75_thr_0.7.mat
xDAEBA = xs;
RDAEBA = R;
load IoU-val-SE_TOrient_YL_Detection_rate_a_0.65_b_0.75.mat
RDREBA = R;
xDREBA=xs;

figure(1); h1 = plot(xDAEB,RDAEB,'g');hold on;box on;grid on;
h2 = plot(xDAYL,RDAYL,'r');
h3 = plot(xDAEBA,RDAEBA,'B');
legend([h1,h2,h3],'Edge Boxes','Replaced Edge Linking','Edge Boxes removed first stage','Location','SouthEast')
xlabel('# of proposals')
ylabel('Detection Rate')
hold off

figure(1); h1 = plot(xDREB,RDREB,'g');hold on;box on;grid on;
h2 = plot(xDRYL,RDRYL,'r');
h3 = plot(xDREBA,RDREBA,'B');
legend([h1,h2,h3],'Edge Boxes','Replaced Edge Linking','Edge Boxes removed first stage','Location','NorthEast')
xlabel('IoU')
ylabel('Detection Rate')
