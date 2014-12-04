clear all;clc;clf
% Using makecircle on all the cov data
yearstart = 2008;
yearend = 2012;
int = '30days_int_';
norm_str = 'norm_cov_data_'; % 'norm...' or 'unnorm...'
N = 7; % Number of 'weeks'
thres1 = 0.975;
thres2 = 0.2;

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
        results = makecircle(cov_data, ii, jj, thres1, thres2, norm);
    end
end