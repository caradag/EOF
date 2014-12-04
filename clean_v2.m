function [station_clean, loggers_clean, x_clean] = clean_v2(station, logger, x)
% v2 has an extra input logger
% Remove the rows and columns in the covariance matrix x of size mxm that contains 
% all NaN's and the corresponding boreholes and loggers that has faulty data. 
% Inputs: station - cell array of size mx1. 
%         logger - cell array of size mx1. 
%         x - covariance matrix of size mxm. 
% Outputs: station_clean - cell array with corresponding boreholes (rows) removed
%          loggers_clean - cell array with corresponding loggers (rows) removed
%          x_clean - covariance matrix with rows and/or columns containing 
%                    all NaN's removed

x_clean = x;
station_clean = station;
loggers_clean = logger;

mat = isnan(x);
vec2 = all(mat,2)==1;  % indices of rows that are all NaN's
vec3 = all(mat,1)==1;  % indices of columns that are all NaN's

x_clean(vec2,:) = [];  % deleting rows
x_clean(:,vec3) = [];  % deleting columns
station_clean(vec3) = [];
loggers_clean(vec3) = [];

end