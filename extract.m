function [p_out, t_out] = extract(tstart, tfinal, t_in, p_in, dt,...
    normalize, max_nan, max_succ_nan, interp_meth)
% Extract the segment of data from column vectors t (time) and p (pressure), 
% remove the NaN's, and interpolating them (Units of all time-related inputs
% and outputs are in days)
% Inputs: tstart, tfinal - scalars in days
%         t_in, p_in - column vectors of equal size   
%         dt - scalar in days, user-specified; default is 1 minute if not 
%              specified
%         normalize - Boolean, user-specified. True normalizes p_out; 
%                     default is false (not normalized)
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
    dt = 1/60/24; 
end

if nargin < 6 || isempty(normalize);
    normalize = false;
end

if nargin < 7 || isempty(max_nan);
    max_nan = 10; 
end

if nargin < 8 || isempty(max_succ_nan);
    max_succ_nan = 3; 
end

if nargin < 9 || isempty(interp_meth);
    interp_meth = 'spline';
end

% Finding the indices of the data points of 1 data point before tstart and 
% 1 data point after tfinal (for interpolation of column vector p)
t1 = t_in(t_in>=tstart); t1 = t1(t1<=tfinal);
n1 = find(t_in==t1(1));    % n1 is the index of the data point of 1 data 
                           % point before tstart in t_in
n2 = find(t_in==t1(end));  % n2 is the index of the data point of 1 data
                           % point after tfinal in t_in
if n1>1;
    n1 = n1-1;
else
end

if n2<length(t_in);
    n2 = n2+1;
else
end

% Extracting the column vectors p and t with 1 extra data point on each ends
p1 = p_in(n1:n2); 
t1 = t_in(n1:n2);

% Removing the NaN's in t1 and p1
nan_num = sum(isnan(p1));
nanspacing = diff(find(isnan(p1)));
if nan_num <= max_nan && sum(nanspacing <= max_succ_nan)<1;
    z = isnan(p1);
    t1(z) = []; p1(z) = [];
end

% Interpolating p onto t with abscissa dt
t_out = (tstart:dt:tfinal).';   
p_out = interp1(t1,p1,t_out,interp_meth);

% Removing the mean and normalizing the column vector p_out
%p_out = p_out - mean(p_out);

if normalize;
    norm_fac = (1/(tfinal-tstart))*(p_out.'*p_out)*dt;
    p_out = p_out./norm_fac;
end
end