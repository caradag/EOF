function [cov, p_mean, C] = cov_v3(p, t, normalize)
%% v3 computes either the normalized or unnormalized covariance matrix  
%  and assumes delta t throughout t is constant
% Inputs: p - matrix of time series (mxn). Each m row represents a time
%             series. Each time series has n points
%         t - column vector of time stamps (nx1)
%         normalize - Boolean, user-specified. True normalizes cov; 
%                     default is false (unnormalized)
% Outputs: cov - square matrix (mxm). It's entries are the corresponding
%                covariances of each time series with each other. 
%          p_mean - column vector (mx1) of the mean of all the time series
%          C - column vector (mx1) of the normalizing factor of all the 
%              time series

% Setting user-specified input to default settings if unspecified
if nargin < 3 || isempty(normalize);
    normalize = false;    
end

% Calculations
A = size(p);  
m = A(1); n = A(2);
p_mean = mean(p,2);
p2 = p - repmat(p_mean,1,n);
pre_cov = p2*(p2.');
C = sqrt((t(2)-t(1))/(t(end)-t(1))*diag(pre_cov));

if normalize;
    C_mat1 = C*C.';    
else
    C_mat1 = ones(m,m);
end

cov = (t(2)-t(1))./(t(end)-t(1))./C_mat1.*pre_cov;

end