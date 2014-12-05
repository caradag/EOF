clear all;clc;clf
% Using makecircle_v4 on all the cov data
tic;
yearstart = 2008; 
yearend = 2014;

startWeek=1;
endWeek=73;

startDOW=1;
endDOW=5; % Should be equal than windowLength if all data is to be processed

yearstart = 2014; 
startWeek=29;
startDOW=1;
yearend = yearstart;
endWeek=startWeek;
endDOW=startDOW;


windowLength=5; % To process data of 5 days windows
windowOffsets=5; % Number of offset in the window for example:
                 %      5 => 5 days window with 1 day time steps
                 %      3 => 15 days window with 5 day time steps
normData = false; % true to process normalized data, false for unnormalized

% Looping through the years
for year=yearstart:yearend
    % Looping through the 'weeks'
    fprintf('Year %d\n',year);
    for week=startWeek:endWeek
        fprintf('   Window %d out of %d\n',week,endWeek-startWeek+1);
        %looping trough folders identified by Day Of Window number
        for dow=startDOW:endDOW
            if normData
                folderName=sprintf('Results/Normalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
                cov_data = sprintf('%ddays_int_norm_cov_data_%d_%d.mat',windowLength,year,week);
            else
                folderName=sprintf('Results/Unnormalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
                cov_data = sprintf('%ddays_int_unnorm_cov_data_%d_%d.mat',windowLength,year,week);
            end
            if exist([folderName cov_data],'file')
                results = makecircle_v4(cov_data, year, week, normData, folderName,'map','off');
                fprintf('      Day of window %d out of %d DONE\n',dow,windowOffsets);
            else
                fprintf('      Day of window %d out of %d DATA FILE DOESN''T EXIST: %s%s\n',dow,windowOffsets,folderName,cov_data);
            end
        end
    end
end
fprintf('DONE in %.1f minutes\n',toc/60);