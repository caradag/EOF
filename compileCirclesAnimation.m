clear all;clc;clf

aviobj=avifile('dancing_circles_full_frame.avi','fps',7,'compression','none');

yearstart = 2008; 
yearend = 2014;

startWeek=1;
endWeek=73;

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
        fprintf('   Window %d out of %d\n',week,numberOfWeeks);
        %looping trough folders identified by Day Of Window number
        for dow=1:windowOffsets
            if normData
                folderName=sprintf('Results/Normalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
            else
                folderName=sprintf('Results/Unnormalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
            end
            filename=sprintf('%d_%d_Map 1.png',year,week);
            fileFullPathName=[folderName 'Plots/' filename];
            if exist(fileFullPathName,'file')
                frame=imread(fileFullPathName);
                %frame=imcrop(frame,[156 67 541 735]);
                aviobj = addframe(aviobj,frame);
            else
                fprintf(2,'File not found : %s\n',fileFullPathName)
            end
        end
    end
end
clo=close(aviobj)