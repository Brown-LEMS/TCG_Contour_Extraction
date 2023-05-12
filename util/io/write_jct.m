function write_jct(filename, G, total_num_edges)
    T_jcts_out = locate_T_jct_from_graph(G, total_num_edges);
    % output conners and junctions
    line_nbr_th = 5;
    fprintf('saving into .jct file\n');
    fout = fopen(filename, 'w');
    fprintf(fout, ' Y-junctions and corners: edge_id [connected contour indices] [assocated contour outward directions]\n\n');
    % locate each junction from graph, compute the direction of each
    % connected contour in order
    for v = 1:length(G.var)
        jct_id = G.var(v).nbrs_fac;
        if(length(jct_id) < 2)
            continue;
        end
        
        jct_contour_dir = zeros(size(jct_id));
        for c = 1:length(jct_id)
            cur_c_id = jct_id(c);
            cur_cfrag = G.fac(cur_c_id).contour;
            % consider to compute contour directions from both ends
            edge_0 = cur_cfrag(1,:);
            edge_1 = cur_cfrag(min(line_nbr_th, size(cur_cfrag, 1)), :);
            if(G.fac(cur_c_id).nbrs_var(2) == G.var(v).id)
                edge_0 = cur_cfrag(end, :);
                edge_1 = cur_cfrag(max(1, size(cur_cfrag, 1) - line_nbr_th), :);
            end
            dx = edge_1(1) - edge_0(1);
            dy = edge_1(2) - edge_0(1);
            jct_contour_dir(c) = atan2(dy, dx);
        end
        % jucntion format: edge_id [ids of connected contours] [direction of connected contours]
        s = [  num2str(G.var(v).actual_edge_id' - 1) ' ' mat2str(jct_id' - 1) ' ' mat2str(jct_contour_dir')];
        fprintf(fout, '%s\n', s);
    end
    
    fprintf(fout, '\n T-junctions: [edge_id stem_contour_id reach_contour_id]\n\n');
    for i = 1:size(T_jcts_out, 1)
        s = mat2str(T_jcts_out(i,:)-1);
        fprintf(fout, '%s\n', s);
    end
    fclose(fout);

end