function results = cleanup(rem_sensors, covdataname, yr, wk, thres1, thres2, norm)
% Output a new cov_data without the information from sensor 
% Inputs: rem_sensor - cell array, the list of sensor information that are going to 
%                      removed.
%         covdataname - string, the name of the cov_data file.
%         yr, wk, thres1, thres2, norm - the inputs in makecircle_v2.
% Output: results - string, irrelevant. See saved cov_data.

load(covdataname);
for ii=1:length(rem_sensors);
    n = find(strcmp(rem_sensors{ii},sensors_clean));
    sensors_clean(n)=[];
    cov_clean(n,:)=[]; cov_clean(:,n)=[];
end
filename = [covdataname ' (2)'];
save(filename, 'cov_clean', 'sensors_clean', 'end_time', 'start_time', 'press_mean');

% Using makecircle_v2 on the new cov data
results = makecircle_v2(cov_data, yr, wk, thres1, thres2, norm);
end