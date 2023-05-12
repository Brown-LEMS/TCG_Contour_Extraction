function [x_coords, y_coords] = interpolate_cfrag(cfrag)

    c_len = contour_length_mex(cfrag');
    pathXY = cfrag(:, 1:2);
    pathXY = unique(pathXY, 'rows', 'stable');
    if(size(pathXY, 1)<2)
        x_coords = pathXY(1);
        y_coords = pathXY(2);
        return;
    end
    
    % interpolate of each contour, approximately at 0.5 step length
    % this fill the gaps
    stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
    stepLengths = [0; stepLengths]; % add the starting point
    cumulativeLen = cumsum(stepLengths);
    finalStepLocs = linspace(0,cumulativeLen(end), round(c_len)*2); 
    finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
    x_coords = finalPathXY(:,1);
    y_coords = finalPathXY(:,2);
%     x_coords = round(finalPathXY(:,1))+1;
%     y_coords = round(finalPathXY(:,2))+1;
%     x_coords = min(x_coords, w);
%     x_coords = max(x_coords, 1);
%     y_coords = min(y_coords, h);
%     y_coords = max(y_coords, 1);