function cov = cov_v1 (x, t)
%% v1 assumes delta t throughout t is constant
% Inputs: x - matrix of time series of size m x n. Each row represents a 
%             time series. Each time series has n points
%         t - column vector of size n x 1
% Output: cov - square matrix of size m x m. It's entries are the
%               corresponding covariances of each time series with each other 

% Checks
A = size(x); % A is row vector [m, n]
n = length(t);

if A(2) ~= n;
    error ('number of time steps in x and t are not the same')
end
    
% Calculations
xmean = mean(x,2);
cov = (x*(x.'))/n - xmean*(xmean.');
end
