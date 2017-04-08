function  [ N ]  = bootstrap_brain_func(F, nummales)
    SET1 = F(1:nummales);
    SET2 = F((nummales + 1):end);
    %Find the mean of the first group and the second group.
    MeanM = mean(SET1);
    MeanF = mean(SET2);
    %Calculate the test statistic
    N = MeanM - MeanF;
end
