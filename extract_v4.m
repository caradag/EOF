function [p_out, t_out] = extract_v4(tstart, tfinal, t_in, p_in, dt,...
    max_nan, max_succ_nan, interp_meth)
%% v4 does not remove the mean and then normalize, and does not compute the
%  mean and normalizing factor
% Extract the segment of data from column vectors t(time) and p(pressure), 
% remove the NaN's, and interpolate them (Units of all time-related 
% inputs and outputs are in days). t_out is made monotonic by removing all 
% "older" overlapping data points and the corresponding data points in 
% p_out is removed as well before being passed for interpolation. 
% Inputs: tstart, tfinal - scalars in days. tstart must NOT be t_in(1) and 
%                          tfinal must NOT be t-in(end)
%         t_in, p_in - column vectors of equal size   
%         dt - scalar in days, user-specified; default is 2 minutes if not 
%              specified
%         max_nan - maximum number of NaN's tolerated in vector p for
%                   removal; default is 10 if not specified 
%         max_succ_nan - maximum number of successive NaN's tolerated in 
%                        vector p for removal; default is 3 if not
%                        specified
%         interp_meth - interpolation method (string); default is 'spline'
%                       if not specified
% Outputs: p_out, t_out - column vectors of equal size with data between 
%                         tstart and tfinal

% Setting user-specified inputs to default if not specified
if nargin < 5 || isempty(dt);
    dt = 2/60/24; 
end

if nargin < 6 || isempty(max_nan);
    max_nan = 10;
end

if nargin < 7 || isempty(max_succ_nan);
    max_succ_nan = 3;
end

if nargin < 8 || isempty(interp_meth);
    interp_meth = 'spline';
end

% Deleting the NaN's in the time stamps and the corresponding location in
% the pressure series
p_in(isnan(t_in)) = []; t_in(isnan(t_in)) = [];

% Finding the indices of the nearest data points before tstart (n_minus) 
% and after tfinal (n_plus) which are not NaN's, for interpolation of  
% column vector p 
n_minus = find((t_in<tstart&~isnan(p_in)), 1, 'last');
n_plus  = find((t_in>tfinal&~isnan(p_in)), 1);

n_inside = find(tstart <= t_in & t_in <= tfinal); % In case n_minus or 
                                                  % n_plus is/are empty

% Extracting the column vectors p and t with 1 extra data point on each 
% ends
n1 = [n_minus; n_inside; n_plus];
t1 = t_in(n1);
p1 = p_in(n1);

% Removing the NaN's in t1 and p1
nan_num = sum(isnan(p1));
nanspacing = diff(find(isnan(p1)));
if nan_num <= max_nan && sum(nanspacing <= max_succ_nan)<1;
    n2 = isnan(p1); % n2 are the locations of the NaN's
    t1(n2) = []; p1(n2) = [];
end


% Making t1 monotonic by removing the "older" overlapping data in both t1
% and p1
[t1,I] = sort(t1);
n3 = find(diff(I)<0)+1; % Indices of the data that need to be removed 
t1(n3) = []; p1(n3) = [];
% Removing the "older" but same time data in both t1 and p1
n4 = find(diff(t1)==0); % Indices of the data that need to be removed
t1(n4) = []; p1(n4) = [];

% Interpolating p onto t with abscissa dt
t_out = (tstart:dt:tfinal).';   

if nan_num > max_nan || sum(nanspacing <= max_succ_nan)>=1 || length(p1)<=1;
    p_out = nan(length(t_out),1); % Setting p_out to be a column vector of 
                                  % NaN's if tolerance(s) is/are exceeded
else
    p_out = interp1(t1,p1,t_out,interp_meth);
end
end