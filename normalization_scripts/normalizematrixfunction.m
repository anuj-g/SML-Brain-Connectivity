function [ N ] = normalizematrixfunction(M)
%Normalize the row by dividing each element by row sum
R = diag(1./sum(M,2))*M;
%Find the maximum weight in this matrix
maxe = max(R(:));
N = R./maxe;
end

