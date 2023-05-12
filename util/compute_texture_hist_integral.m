function [texton_hist_left_vector, texton_hist_right_vector, texton_hist_left_integral, texton_hist_right_integral] = compute_texture_hist_integral(x,y,N, nbr_width, tmap)

    % lateral edge sparsity
    [h,w, ~] = size(tmap);
    
    % N's direction matters
    left_x = round(x -  nbr_width*N(:,1));
    left_y = round(y -  nbr_width*N(:,2));
    right_x = round(x +  nbr_width*N(:,1));
    right_y = round(y +  nbr_width*N(:,2));

    % deal with the local exceeding the border
    left_x = max(left_x, ones(size(left_x)));
    left_y = max(left_y, ones(size(left_y)));
    left_x = min(left_x, w*ones(size(left_x)));
    left_y = min(left_y, h*ones(size(left_y)));
    right_x = max(right_x, ones(size(right_x)));
    right_y = max(right_y, ones(size(right_y)));
    right_x = min(right_x, w*ones(size(right_x)));
    right_y = min(right_y, h*ones(size(right_y)));
    
    texton_hist_left_vector = zeros(64, length(x));
    texton_hist_right_vector = zeros(64, length(x));
    
    for i = 1:length(x)
        x_coods_left = linspace(left_x(i), x(i), nbr_width+1);
        y_coods_left = linspace(left_y(i), y(i), nbr_width+1);
        lin_ind_left = sub2ind([h,w], round(y_coods_left(1:end-1)), round(x_coods_left(1:end-1)));
        texton_hist_left_vector(:,i) = (hist(tmap(lin_ind_left), 1:64))';

        x_coods_right = linspace(right_x(i), x(i), nbr_width+1);
        y_coods_right = linspace(right_y(i), y(i), nbr_width+1);
        lin_ind_right = sub2ind([h,w], round(y_coods_right(1:end-1)), round(x_coods_right(1:end-1)));
        texton_hist_right_vector(:,i) = (hist(tmap(lin_ind_right), 1:64))';
        
        
    end
    
    texton_hist_left_integral = cumsum(texton_hist_left_vector, 2);
    texton_hist_right_integral = cumsum(texton_hist_right_vector, 2);
    

end