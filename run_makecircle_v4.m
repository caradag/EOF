clear all;clc;clf
% Using makecircle_v4 on all the cov data
yearstart = 2008; 
yearend = 2014;
int = '5days_int_';
norm_str = 'unnorm_cov_data_'; % 'norm...' or 'unnorm...'
N = 73; % Number of 'weeks'

if strcmp(norm_str(1),'n');
    norm = true; 
else
    norm = false;
end

% Looping through the years
for ii=yearstart:yearend;
    % Looping through the 'weeks'
    for jj=1:N;
        cov_data = [int norm_str num2str(ii) '_' num2str(jj) '.mat'];
        results = makecircle_v4(cov_data, ii, jj, norm);
    end
end