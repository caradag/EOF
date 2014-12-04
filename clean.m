function [station_clean, x_clean] = clean(station, x)
% Remove the rows and columns in the covariance matrix x of size mxm that 
% contains all NaN's and the corresponding boreholes that has faulty data. 
% Inputs: station - cell array or column vector of size mx1. 
%         x - covariance matrix of size mxm. 
% Outputs: station_clean - cell array or column vector with corresponding 
%                          boreholes (rows) removed
%          x_clean - covariance matrix with rows and/or columns containing 
%                    all NaN's removed

x_clean = x;
station_clean = station;

mat = isnan(x);
vec2 = all(mat,2)==1;  % indices of rows that are all NaN's
vec3 = all(mat,1)==1;  % indices of columns that are all NaN's

x_clean(vec2,:) = [];  % deleting rows
x_clean(:,vec3) = [];  % deleting columns
station_clean(vec3) = [];

end