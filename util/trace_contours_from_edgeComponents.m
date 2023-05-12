function [edginfo, edgelist2, cfrags_idx] = trace_contours_from_edgeComponents (edgeComponents, E, O)
    edgelist2 = cell(1,0);
    cfrags_idx = cell(1,0);
    [yy,xx] = find(edgeComponents~=0);
    [h, w] = size(E);
    Angle = O(sub2ind([h,w], yy, xx));
    mag = E(sub2ind([h,w], yy, xx));
    edginfo  = [xx yy Angle mag];
    % construct edge idx map
    idx_map = zeros(size(E));
%     idx_map(sub2ind([h,w], edginfo(:,2), edginfo(:,1))) = 1:size(edginfo, 1);
    for eid = 1:size(edginfo, 1)
        idx_map(edginfo(eid, 2), edginfo(eid, 1)) = eid;
    end
    
    component_ids = unique(edgeComponents);
    cnt = 0;
%     figure(2); hold on;
    for i = 1:length(component_ids)
       
        id = component_ids(i);
        if(id ==0)
            continue;
        end
        
        % extracted a connected component -> it missed the shared node
        [e_Y, e_X] = find(edgeComponents==id);
        % remove the extremely small group
        if(length(e_Y)<=3)
            edgeComponents(sub2ind([h,w], e_Y, e_X)) = 0;
            continue;
        end
       
        % local trace
        edges = [];
        dx = [-1, 0, 1,-1, 1,-1, 0, 1];
        dy = [-1,-1,-1, 0, 0, 1, 1, 1];
        
        start = [e_X(1), e_Y(1)];
        edges = [edges; [start O(start(2), start(1)) E(start(2), start(1))]];
        edgeComponents(start(2),start(1)) = 0;
        
        nbr = find(edgeComponents(sub2ind([h,w], dy+start(2),dx+start(1)))==id);
        dist = dy(nbr).^2 + dx(nbr).^2;
        [~, sort_idx] = sort(dist);
        nbr = nbr(sort_idx);

        nbr1 = nbr(1);
        prev1 = start;
        nbr2 = [];
        if(length(nbr)>1) % it's possible start is an end
            nbr2 = nbr(2);
        end
        prev2 = start;
        
        while(~isempty(nbr1))
            % if there are two neighbors use the closer one
            if(length(nbr1)>1)
                dist = dy(nbr1).^2 + dx(nbr1).^2;
                [~, min_idx] = min(dist);
                nbr1 = nbr1(min_idx);
            end
            
            cur1 = prev1 + [dx(nbr1) dy(nbr1)];
            if(edgeComponents(cur1(2),cur1(1)) == 0)
                break;
            end
            edges = [edges; [cur1 O(cur1(2),cur1(1)) E(cur1(2),cur1(1))]];
            edgeComponents(cur1(2),cur1(1)) = 0;
            prev1 = cur1;
            nbr1 = find(edgeComponents(sub2ind([h,w], dy+prev1(2), dx+prev1(1)))==id);
        end
        
        while(~isempty(nbr2))
            % if there are two neighbors use the closer one
            if(length(nbr2)>1)
                dist = dy(nbr2).^2 + dx(nbr2).^2;
                [~, min_idx] = min(dist);
                nbr2 = nbr2(min_idx);
            end
            cur2 = prev2 + [dx(nbr2) dy(nbr2)];
            if(edgeComponents(cur2(2),cur2(1)) == 0)
                break;
            end
            edges = [[cur2 O(cur2(2),cur2(1)) E(cur2(2),cur2(1))]; edges];
            edgeComponents(cur2(2),cur2(1)) = 0;
            prev2 = cur2;
            nbr2 = find(edgeComponents(sub2ind([h,w], dy+prev2(2),dx+prev2(1)))==id);
        end
        
        cnt =  cnt+1;
%         plot(edges(:,1), edges(:,2), '-o');
        edgelist2{cnt} = edges;
        cfrags_idx{cnt} = idx_map(sub2ind([h,w], edges(:,2), edges(:,1)))';
    end
    
end
