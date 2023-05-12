function LabelEdgemap = convert_cfrags_to_EdgeGroupMap (cfrags, h, w, interp)
if(nargin<4)
    interp = 1;
end
ratio = 2;

LabelEdgemap = zeros(h,w);
for i = 1:length(cfrags)

    c_len = contour_length_mex(cfrags{i}');
    pathXY = cfrags{i}(:, 1:2);
    pathXY = unique(pathXY, 'rows', 'stable');
    if(size(pathXY, 1)<2)
        continue;
    end
    if(~interp)
        x_coords = round(pathXY(:,1))+1;
        y_coords = round(pathXY(:,2))+1;
        x_coords = min(x_coords, w);
        x_coords = max(x_coords, 1);
        y_coords = min(y_coords, h);
        y_coords = max(y_coords, 1);

        LabelEdgemap(sub2ind([h,w], y_coords, x_coords)) = i;
    else
        % interpolate of each contour, approximately at 0.5 step length
        % this fill the gaps
        stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
        stepLengths = [0; stepLengths]; % add the starting point
        cumulativeLen = cumsum(stepLengths);
        finalStepLocs = linspace(0,cumulativeLen(end), round(c_len)*ratio); 
        finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
        x_coords = round(finalPathXY(:,1))+1;
        y_coords = round(finalPathXY(:,2))+1;
        x_coords = min(x_coords, w);
        x_coords = max(x_coords, 1);
        y_coords = min(y_coords, h);
        y_coords = max(y_coords, 1);
%         LabelEdgemap(sub2ind([h,w], y_coords, x_coords)) = i;

        binary_map = zeros(h,w);
        binary_map(sub2ind([h,w], y_coords, x_coords)) = i;
%         binary_map = bwmorph(binary_map,'clean');     % Remove isolated pixels
        binary_map = bwmorph(binary_map,'skel',Inf);  % and make sure edges are thinned. I

        LabelEdgemap(binary_map>0) = i;

    end
end

end