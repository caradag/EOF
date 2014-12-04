clear all;clc
% Creating 'clean data v6.mat' (To update sensors and sensor_count, change S11P21 
% grid to 'R20C21.25') 
load('clean data v5.mat');
sensors = fieldnames(data);
sensor_count = length(sensors);
data.S11P21.grid={'R20C21.25'};
filename = 'clean data v6';
save(filename,'data','sensors','sensor_count');