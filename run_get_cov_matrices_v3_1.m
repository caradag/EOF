clear all;clc
% Using get_cov_matrices to save the covariance matrices
datatype='pressure';
yearstart='2008';
yearend='2014';
tstartinyear='5-Jan-2000';
tendinyear='31-Dec-2000';
t_interval=5;
folder = 'Unnormalized press, 5 days int, 2008 to 2014 (5)';
dt=2/60/24;
max_nan=10;
max_succ_nan=4;
interp_meth='spline';
max_dt=20/60/24;
max_miss=max_succ_nan;
max_miss_int=max_succ_nan;
max_miss_out=max_succ_nan;
normalize=false;
results = get_cov_matrices_v3_1(datatype, yearstart, yearend,...
    tstartinyear, tendinyear, t_interval, folder, dt, max_nan,...
    max_succ_nan, interp_meth, max_dt, max_miss, max_miss_int,...
    max_miss_out, normalize);