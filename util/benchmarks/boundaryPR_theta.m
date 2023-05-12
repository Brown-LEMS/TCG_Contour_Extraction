function [thresh,cntR,sumR,cntP,sumP] = boundaryPR_theta(pb,pb_orient,segs,orient,theta_diff)
% function [thresh,cntR,sumR,cntP,sumP] = boundaryPR(pb,segs,nthresh)
%
% Calcualte precision/recall curve.
% If pb is binary, then a single point is computed.
% The pb image can be smaller than the segmentations.
%
% INPUT
%	pb		Soft or hard boundary map.
%	segs		Array of segmentations.
%	[nthresh]	Number of points in PR curve.
%
% OUTPUT
%	thresh		Vector of threshold values.
%	cntR,sumR	Ratio gives recall.
%	cntP,sumP	Ratio gives precision.
%
% See also boundaryPRfast.
% 
% David Martin <dmartin@eecs.berkeley.edu>
% January 2003

bmap = double(bwmorph(pb,'thin',inf));

% Numb ground truth segmentations
nsegs = length(segs);

% zero all counts
cntR = 0;
sumR = 0;
cntP = 0;
sumP = 0;
[height,width] = size(pb);
% compute boundary maps from segs
bmaps = cell(size(segs));
for i = 1:nsegs,
  bmaps{i} = double(seg2bmap(segs{i},width,height));
end

% make sure the boundary maps are thinned to a standard thickness
for i = 1:nsegs,
  bmaps{i} = bmaps{i} .* bwmorph(bmaps{i},'thin',inf);
end

% accumulate machine matches, since the machine pixels are
% allowed to match with any segmentation
accP = zeros(size(pb));
% compare to each seg in turn

for i = 1:nsegs
    
    % compute the correspondence
    [match1,match2] = correspondPixels_theta(bmap,bmaps{i},pb_orient,orient{i},theta_diff);
    %[match1,match2] = correspondPixels(bmap,bmaps{i});
    
    % accumulate machine matches
    accP = accP | match1;
    
    % compute recall
    sumR = sumR + sum(bmaps{i}(:));
    cntR = cntR + sum(match2(:)>0);
    
end
% compute precision
sumP = sumP + sum(bmap(:));
cntP = cntP + sum(accP(:));

thresh=1;