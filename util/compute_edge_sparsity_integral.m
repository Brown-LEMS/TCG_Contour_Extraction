function [edge_sparcity_vector, edge_sparsity_integral] = compute_edge_sparsity_integral(x,y,N, nbr_width, edge_map)
    % lateral edge sparsity
    [h,w, ~] = size(edge_map);
    edge_map = edge_map>0;
    
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
    
    edge_sparcity_vector = zeros(size(x));
    
    for i = 1:length(x)
        x_coods = linspace(left_x(i), right_x(i), 2*nbr_width+1);
        y_coods = linspace(left_y(i), right_y(i), 2*nbr_width+1);
        lin_ind = sub2ind([h,w], round(y_coods), round(x_coods));
        edge_sparcity_vector(i) = sum(edge_map(lin_ind));
        
    end
    
    edge_sparsity_integral = cumsum(edge_sparcity_vector);    
end