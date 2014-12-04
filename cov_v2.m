function cov = cov_v2 (x, t, tol_mean)
%% v2 defines cov(x_i,x_j)= 1/(t_final - t_initial) int_(t_initial)^(t_final) x_i(t) x_j(t) dt 
%  and assumes (i) all time series in matrix x have approximately zero mean (see tol_mean), 
%          and (ii) delta t throughout t is constant
% Inputs: x - matrix of time series of size m x n. Each row represents a 
%             time series. Each time series has n points
%         t - column vector of size n x 1
%         tol_mean - tolerated error in each time series means from zero. If tol_mean is not 
%                    specified, default value is sqrt(eps)
% Output: cov - square matrix of size m x m. It's entries are the
%               corresponding covariances of each time series with each other 

% Checks
A = size(x); % A is row vector [m, n]
n = length(t);
xmean = mean(x,2);

if nargin < 3
    tol_mean = sqrt(eps);
end

if sum(xmean > tol_mean*ones(A(1),1))~=0;
    warning ('Some or all the means of the time series in matrix x are not within the tolerated error from zero')
end
    
if A(2) ~= n;
    error ('Number of time steps in x and t are not the same')
end

% Calculations
cov = ((t(2)-t(1))/(t(end)-t(1)))*x*(x.');
end