clear all;clc
% Using get_cov_matrices to save the covariance matrices into the 
% CURRENT FOLDER
datatype='temperature';
yearstart='2008';
yearend='2014';
tstartinyear='1-Jan-2000';
tendinyear='31-Dec-2000';
t_interval=10;
dt=2/60/24;
max_nan=10;
max_succ_nan=4;
interp_meth='spline';
max_dt=20/60/24;
max_miss=max_succ_nan;
max_miss_int=max_succ_nan;
max_miss_out=max_succ_nan;
normalize=true;
results = get_cov_matrices_v3(datatype, yearstart, yearend,...
    tstartinyear, tendinyear, t_interval, dt, max_nan, max_succ_nan,...
    interp_meth, max_dt, max_miss, max_miss_int, max_miss_out, normalize);