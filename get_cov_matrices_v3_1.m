function results = get_cov_matrices_v3_1(datatype, yearstart, yearend,...
    tstartinyear, tendinyear, t_interval, folder, dt, max_nan,...
    max_succ_nan, interp_meth, max_dt, max_miss, max_miss_int,...
    max_miss_out, normalize)
% Operates on the data within yearstart and yearend. Covariance matrices
% are calculated over the time interval, t_interval from tstartinyear till 
% tendinyear. The clean covariance matrix (cov_clean), means (press_mean),
% normalizing factors (norm_fac), clean sensors (sensors_clean), the
% starting and ending time (start_time, end_time) over which the 
% covariance matrix is calculated are SAVED into the CURRENT FOLDER. 
% [Note: Units of t_interval and dt are in days.]
% v3_1 has an extra input 'folder'
%
% Inputs: datatype - string, 'pressure' or 'temperature' ONLY.
%         yearstart, - strings, the starting and ending year over which 
%         final_time   the data are considered. Eg: '2008'
%         tstartinyear, - strings, the start and end time within the year.
%         tendinyear      Eg: '30-May-2000' to '31-Dec-2000' (The year
%                         doesn't matter but both must be the same).
%                         NOTE: There must be EXACTLY n*t_interval days 
%                         between tstartinyear and tendinyear where n is an
%                         integer
%         t_interval - scalar, the time interval over which each covariance 
%                      matrix is being calculated
%         folder - string, the directory to which the files are saved to.
%         (For inputs dt to max_miss_out, see extract_v5)
%         (For input normalize, see cov_v3)
% Output: results - string, irrelevant. See saved files.

% Setting user-specified inputs to default if not specified
if nargin < 8 || isempty(dt);
    dt = 2/60/24; 
end

if nargin < 9 || isempty(max_nan);
    max_nan = 10; 
end

if nargin < 10 || isempty(max_succ_nan);
    max_succ_nan = 3; 
end

if nargin < 11 || isempty(interp_meth);
    interp_meth = 'spline';
end

if nargin < 12 || isempty(normalize);
    normalize = false;
end

data = load('data 2014 v5 good only.mat');

% Setting up some useful variables
sensors = fieldnames(data);
sensor_count = length(sensors);
startmonth = datestr(tstartinyear,3); % String, starting month in a year
startday = datestr(tstartinyear,7); % String, starting day of the month

numdaysinyear = datenum(tendinyear)-datenum(tstartinyear);
N1 = numdaysinyear/t_interval; % Number of t_interval days in a year
N2 = eval(yearend)-eval(yearstart)+1; % Number of N1's in total
year = eval(yearstart);

% Looping through each year from init_time to final_time
for kk=1:N2
    tstart = datenum(strcat(startday,'-',startmonth,'-',num2str(year)));
    % Looping through every t_interval within a year
    for jj=1:N1;
        tfinal = tstart+t_interval;
        
        % Preallocating the matrix for cov_v3
        press_mat = ones(sensor_count,length(tstart:dt:tfinal));
        loggers = cell(sensor_count,1);
       
        % Looping through the boreholes
        for ii=1:sensor_count;
            % Assigning data as inputs for extract_v4
            t_in = eval(strcat('data.',sensors{ii},'.time.serialtime'));
            if strcmp(datatype,'pressure');
                p_in = eval(strcat('data.',sensors{ii},'.pressure'));
            else
                p_in = eval(strcat('data.',sensors{ii},'.temperature'));
            end
            l = eval(strcat('data.',sensors{ii},'.logger'));
            loggers{ii} = l{1};
            p_in = p_in{1};
            tlim = [min(t_in) max(t_in)];
            
            if tstart<tlim(1) || tfinal>tlim(2);
                % Making all the entries in row ii of press_mat as NaN's if
                % tstart<tlim(1) or tfinal>tlim(2)
                t_out = (tstart:dt:tfinal).';
                press_mat(ii,:) = nan(1,length(t_out));
            else
                % Calling extract_v5
                [p_out, t_out] = extract_v5(tstart, tfinal, t_in, p_in,...
                    dt, max_nan, max_succ_nan, interp_meth, max_dt,...
                    max_miss, max_miss_int, max_miss_out);
                press_mat(ii,:) = p_out.';
            end   
        end
        
        % Calling cov_v3
        [cov, press_mean, norm_fac] = cov_v3(press_mat, t_out, normalize); 
        
        % Calling clean
        [sensors_clean, loggers_clean, cov_clean] = clean_v2(sensors.',...
            loggers, cov);
        
        % Converting serial time to human-readable time and saving the 
        % relevant data
        start_time = datestr(tstart);
        end_time = datestr(tfinal);
        
        if normalize;
            filename = strcat(num2str(t_interval),'days_int','_norm_',...
                'cov_data_',num2str(year),'_',num2str(jj));
        else
            filename = strcat(num2str(t_interval),'days_int','_unnorm_',...
                'cov_data_',num2str(year),'_',num2str(jj));
        end
        
        if strcmp(datatype,'pressure');
            save(['Results/' folder '/' filename ],...
                'cov_clean','press_mean','sensors_clean',...
                'loggers_clean','start_time','end_time');
        else
            temp_mean = press_mean;
            save(['Results/' folder '/' filename],...
                'cov_clean','temp_mean','sensors_clean',...
                'loggers_clean','start_time','end_time');
        end
        
        % Resetting tfinal as the new tstart after saving all relevant data
        tstart = tfinal;
    end
    year = year+1;
end

results='Check saved files';

end