function covStack = getCovMatrices(windowSize, timeStep, doNormalize, doDetrend,options)
% Compute covariance matrices over a moving window trough the dataset
%
% Inputs: 
%   windowSize  - Size of the moving window in days
%                 Default: NONE it is mandatory Data type: Double
%   timeStep    - Amount of time in days to shift the window on each step
%                   Default: 1 day  Data type: Double
%   doNormalize - If apply normalization to the data amplitudes or not
%                   Default: false  Data type: boolean
%   doDetrend   - If apply detrending to the data time series or not
%                   Default: false  Data type: boolean
%   options     - Optional structure with computation options.
%                 It can include any of the following fields
%                  
%    datatype  > Either 'pressure' or 'temperature'
%                Default: 'pressure' Data type: char string
%    startTime > Star time of the computation.
%                Default: start of dataset. Data type: matlab serial time
%    endTime   > End time of the computation.
%                Default: end of dataset. Data type: matlab serial time
%    dt        > Time step for data interpolation to a regular sampling
%                Default: 2 minutes Data type: double in DAYS
%  interp_meth > Method used for interpolation as defined for interp1
%                Default: 'linear' Data type: char string
%  save_to_disk> Wether save data to disk (one mat file per window) or not
%                Default: fase if there is output argument in the call
%                         true otherwise. Data type: boolean  
%  folder      > Folder were data files will be saved. If folder doesn't
%                exist it will be created.
%                Default: current folder Data type: char string
%  max_nan     > Option pass to extract_v5. Default 10
%  max_succ_nan> Option pass to extract_v5. Default 4
%  max_dt      > Option pass to extract_v5. Default 20 min
%  max_miss    > Option pass to extract_v5. Default 4
%  max_miss_int> Option pass to extract_v5. Default 4
%  max_miss_out> Option pass to extract_v5. Default 4
%
% Output: If there is an output argument in the cal it will output a vector
%         of sctructures containing the covariances and other relevant data
%         in the following fields:
%  cov         > Covariance matrix
%  timeLims    > Start and end time of the time window (matlab timestamp)
%  sensors     > Sensor IDs of the ones used for computation (cell array)
%  loggers     > Logger ID for each sensor. (cell array)
%  n           > Number of sensors used in computation (integer)
%  window      > Correlative number of the time window (integer)
%
    %% DEALING WITH INPUTS AND ASSINGNING DEFAULTS VALUES IF NEEDED
    if nargin<1
        error('GET_COV_MAT:No_window_length','Window length in days MUST be specified');
    end
    if nargin<2
        timeStep=1;
    end
    if nargin<3
        doNormalize=false;
    end
    if nargin<4
        doDetrend=false;
    end

    if nargin<5
        options=struct;
    end
    if ~isfield(options,'datatype')
        options.datatype='pressure';
    end
    if ~isfield(options,'startTime')
        options.startTime=-Inf; %[2008 1 1];
    end
    if ~isfield(options,'endTime')
        options.endTime=Inf; %[2014 12 31];
    end
    if ~isfield(options,'dt')
        options.dt=2/1440;
    end
    if ~isfield(options,'max_nan')
        options.max_nan=10;
    end
    if ~isfield(options,'max_succ_nan')
        options.max_succ_nan=4;
    end
    if ~isfield(options,'interp_meth')
        options.interp_meth='linear';
    end
    if ~isfield(options,'max_dt')
        options.max_dt=20/1440;
    end
    if ~isfield(options,'max_miss')
        options.max_miss=options.max_succ_nan;
    end
    if ~isfield(options,'max_miss_int')
        options.max_miss_int=options.max_succ_nan;
    end
    if ~isfield(options,'max_miss_out')
        options.max_miss_out=options.max_succ_nan;
    end
    if ~isfield(options,'folder')
        options.folder = pwd;
    else
        if ~exist(options.folder,'dir')
            mkdir(options.folder);
        end
    end
    if ~isfield(options,'save_to_disk')    
        if nargout>0
            options.save_to_disk=false;
        else
            options.save_to_disk=true;
        end
    end
    covStack=[];
    

    %% LOADING DATA AND FINDING DEFAULT TIME LIMITS IF NEEDED
    disp('Loading data file...');
    data = load('data 2014 v5 good only.mat');

    % Setting up some useful variables
    sensors = fieldnames(data);
    sensor_count = length(sensors);
    % If start and end time are not finite we loop trough the data to find the
    % minimum and maximum time stamps
    if ~isfinite(options.startTime) || ~isfinite(options.endTime)
        mint=Inf;
        maxt=-Inf;
        for s=1:sensor_count
            mint=min(mint,data.(sensors{s}).time.serialtime(1));
            maxt=max(maxt,data.(sensors{s}).time.serialtime(end));
        end
        if ~isfinite(options.startTime)
            options.startTime=mint;
        end
        if ~isfinite(options.endTime)
            options.endTime=maxt;
        end
    end
    options.startTime=floor(options.startTime);
    options.endTime=ceil(options.endTime);
    windowSize=round(windowSize);

    %% MAIN LOOP TO COMPUTE COVARIANCE MATRICES
    disp('Computing covariances...');
    tic;
    startTimes=options.startTime:timeStep:(options.endTime-timeStep);
    nWindows=length(startTimes);
    nSamples=floor(windowSize/options.dt)+1;

    % Looping through each time window
    for  windowCount=1:nWindows
        tstart =startTimes(windowCount);
        tfinal = tstart+windowSize;
        
        messageLength=fprintf('%.1f%% %s',100*windowCount/nWindows,datestr(tstart));

        sensorsInRange=false(1,sensor_count);
        % Looping through the boreholes
        for ii=1:sensor_count;
            tlim = data.(sensors{ii}).time.serialtime([1 end]);
            % If sensor data covers the whole window
            if tlim(1)<=tstart && tlim(2)>=tfinal;
                sensorsInRange(ii)=true;
            end   
        end
        sensorsInRange=find(sensorsInRange);
        nSensorsInRange=length(sensorsInRange);
        % Preallocating the matrix for cov_v3
        press_mat = nan(nSamples,nSensorsInRange);
        loggers = cell(nSensorsInRange,1);

        % Looping through the boreholes
        for ii=1:nSensorsInRange;
            sensorIdx=sensorsInRange(ii);
            % Assigning data as inputs for extract_v4
            t_in = data.(sensors{sensorIdx}).time.serialtime;
            switch options.datatype
                case 'pressure'
                    p_in = data.(sensors{sensorIdx}).pressure{1};
                case 'temperature'
                    p_in = data.(sensors{sensorIdx}).temperature{1};
                otherwise
                    error('GET_COV_MAT:Unknown_data_type','Invalid datatype, you can choose pressure or temperature only.');
            end
            loggers{ii} = data.(sensors{sensorIdx}).logger{1};

            % Calling extract_v5
            [p_out, ~] = extract_v5(tstart, tfinal, t_in, p_in,...
                options.dt, options.max_nan, options.max_succ_nan, options.interp_meth, options.max_dt,...
                options.max_miss, options.max_miss_int, options.max_miss_out);
            
            press_mat(:,ii) = p_out(:);
        end
        sensorIDs=sensors(sensorsInRange);

        % Removing nan data series that might be returned if there number
        % and distribution of NaNs go pass specified tresholds.
        nanCols=all(isnan(press_mat));
        if any(nanCols)
            fprintf('%c',8*ones(messageLength,1));  
            fprintf('%s from %s to %s -> Sensors %s eliminated due to NaNs content\n',datestr(tstart,'yyyy'),datestr(tstart,'mmm-dd'),datestr(tfinal,'mmm-dd'),strjoin(sensorIDs(nanCols),','));  
            messageLength=0;
        end
        
        press_mat=press_mat(:,~nanCols);
        sensorIDs=sensorIDs(~nanCols);
        sensorsInRange=sum(~nanCols);
        
        if doDetrend
            press_mat = detrend(press_mat);
        end
        % Computing covarianve        
        covariance=computeCov(press_mat, doNormalize);

        descriptionTxt='';
        if doNormalize
            descriptionTxt='_Normalized';
        end
        if doDetrend
            descriptionTxt=[descriptionTxt '_Detrended'];
        end
        
        % Saving and/or storing data
        if options.save_to_disk
            % Converting serial time to human-readable time vector
            start_time = datevec(tstart);
            end_time = datevec(tfinal);
            filename=sprintf('cov_%s_%04d-%02d-%02d_%dDays%s.mat',options.datatype,start_time,windowSize,descriptionTxt);
            save([options.folder filesep filename ],'covariance','press_mean','sensorIDs','loggers','start_time','end_time','windowCount');
        end
        if nargout>0
            covStack(end+1).cov=covariance;
            covStack(end).timeLims=[tstart tfinal];                    
            covStack(end).sensors=sensorIDs;                    
            covStack(end).loggers=loggers;                    
            covStack(end).n=sensorsInRange;
            covStack(end).window=windowCount;
        end
        fprintf('%c',8*ones(messageLength,1));  
    end
    fprintf('Done in %.1 minutes',toc/60);
end

function covariance=computeCov(x,normalize)
    % Get the number of samples in the time serie in m
    m = size(x,1);
    % Remove the mean of each colum
    x = bsxfun(@minus,x,sum(x,1)/m);
    % Computing unnormalized covariance
    unnormCov = x' * x;

    if normalize
        % This is equivalent to normalize each time series before
        % by dividing it by the square root of the variance
        C = sqrt(diag(unnormCov));
        normDenominator  = C * C';    
    else
        % The time normalization denominator is
        % normDenominator = TimeSpan/TimeStep
        % Which for regular time step is equal to
        normDenominator= m-1;
    end
    covariance = unnormCov./normDenominator;
end