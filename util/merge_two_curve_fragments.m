function c_merge = merge_two_curve_fragments(c1, c2, node)


% make the direction c1-->node-->c2
if(c1(1,1:2)==node(1:2))
    c1 = fliplr(c1');
    c1 = c1';
end

if(c2(end,1:2)==node(1:2))
    c2 = fliplr(c2');
    c2 = c2';
end

c_merge = [c1; c2(2:end, :)];

end