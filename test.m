for i = 1:70
    if(bootstrap_ccf_node('datasets/set1/normalized', i, 0.01) > 0)
        disp(i);
    end
end