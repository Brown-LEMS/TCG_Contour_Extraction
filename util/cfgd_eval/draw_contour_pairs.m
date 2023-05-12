function draw_contour_pairs(contours1, contours2, thresh, random, col)

if (nargin<3), thresh = 0; end
if (nargin<4), random=1; end
if (nargin<5), col = [0 0 0]; end

con_cnt1 = length(contours1);
con_cnt2 = length(contours2);

% this is for multicolor
colourmp1 = hsv(con_cnt1);    % HSV colour map with con_cnt entries
colourmp2 = hsv(con_cnt2); 
% colourmp1 = colourmp1(randperm(con_cnt1),:);  % Random permutation
% colourmp2 = colourmp2(randperm(con_cnt2),:);  % Random permutation

r = rand(1,1);
g = rand(1,1);
b = rand(1,1);
% this is for map1
colourmp1(:,1) = r;
colourmp1(:,2) = g;
colourmp1(:,3) = b;

% % this is for map2
colourmp2(:,1) = r*0.7;
colourmp2(:,2) = g*0.7;
colourmp2(:,3) = b*0.7;

for i = 1:con_cnt1
    if (size(contours1{i},1)<thresh)
        continue;
    end

    if (random==1)
        line(contours1{i}(:,1)+1, contours1{i}(:,2)+1,'color',colourmp1(i,:), 'LineWidth', 2);
    else

        line(contours1{i}(:,1)+1, contours1{i}(:,2)+1,'color', col);

    end
end

for i = 1:con_cnt2
    if (size(contours2{i},1)<thresh)
        continue;
    end

    if (random==1)
        line(contours2{i}(:,1)+1, contours2{i}(:,2)+1,'color',colourmp2(i,:), 'LineWidth', 1);
    else

        line(contours2{i}(:,1)+1, contours2{i}(:,2)+1,'color', col);

    end
end