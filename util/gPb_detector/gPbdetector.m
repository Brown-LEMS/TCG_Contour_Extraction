function [f] = gPbdetector(im)
%% mPb
[mPb, mPb_rsz, bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = multiscalePb(im, 1);

%% sPb
% outFile2 = strcat(outFile, '_pbs.mat');
[sPb] = spectralPb(mPb_rsz);
% delete(outFile2);

bg1 = max(bg1, [], 3);
bg2 = max(bg2, [], 3);
bg3 = max(bg3, [], 3);
cga1 = max(cga1, [], 3);
cga2 = max(cga2, [], 3);
cga3 = max(cga3, [], 3);
cgb1 = max(cgb1, [], 3);
cgb2 = max(cgb2, [], 3);
cgb3 = max(cgb3, [], 3);
tg1 = max(tg1, [], 3);
tg2 = max(tg2, [], 3);
tg3 = max(tg3, [], 3);
sPb = max(sPb, [], 3);

bg1 = bg1(:);
bg2 = bg2(:);
bg3 = bg3(:);
cga1 = cga1(:);
cga2 = cga2(:);
cga3 = cga3(:);
cgb1 = cgb1(:);
cgb2 = cgb2(:);
cgb3 = cgb3(:);
tg1 = tg1(:);
tg2 = tg2(:);
tg3 = tg3(:);
sPb = sPb(:);

f = [ones(size(bg1)), bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, sPb]';

% [cg,tg] = detCGTG(im);
% l = max(squeeze(cg(:,:,1,:)),[],3);
% a = max(squeeze(cg(:,:,2,:)),[],3);
% b = max(squeeze(cg(:,:,3,:)),[],3);
% t = max(tg,[],3);
% l = l(:);
% a = a(:);
% b = b(:);
% t = t(:);
% f = [ ones(size(b)) l a b t ]';

end
