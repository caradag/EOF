clear all;clc
% Using get_cov_matrices to save the covariance matrices into the CURRENT FOLDER
yearstart='2009';
yearend='2009';
tstartinyear='30-Jun-2000';
tendinyear='30-Jul-2000';
t_interval=5;
dt=2/60/24;
max_nan=30;
max_succ_nan=3;
interp_meth='spline';
normalize=true;
results = get_cov_matrices(yearstart, yearend, tstartinyear, tendinyear,...
    t_interval, dt, max_nan, max_succ_nan, interp_meth, normalize)