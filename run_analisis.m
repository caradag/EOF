%Script to run the different rutines of the EOF analisis

dataSetDescription='15Day_normalized_detrended';

% MAIN OPTIONS
% Covariances computation
windowSize=15; % Day ex t_interval
timeStep=5; % Days
doNormalize=true;
doDetrend=true;

% Animation
doCreateCovAnimation=true;
doCreateCirclesAnimation=false;
covThreshold=0.9;

%% CREATING COVARIANCE MATRICES
% Using getCovMatrices to save the covariance matrices

dataFolder=['Results' filesep dataSetDescription];
covDataFile=['covariances_' dataSetDescription '.mat'];

if exist([dataFolder filesep covDataFile],'file')
    load([dataFolder filesep covDataFile]);
else
    disp(['Covariances data file: ' dataFolder filesep covDataFile ' NOT FOUND']);
    disp('Computing covariances...');
    
    % ADDITIONAL OPTIONS
    options=struct;
    options.save_to_disk=false;
    options.folder = 'test';

    covStack = getCovMatrices(windowSize, timeStep, doNormalize, doDetrend,options);

    % SAVING DATA
    % Optiaonally we can save tha covariance data for later use
    save([dataFolder filesep covDataFile],'covStack');
end
%% CREATING EOF CIRCLE PLOTS
% makeCircles will create a figure for EOFs as a circle plot, one for each
% time window and each one of the the mos significant eigenvectors
% Figures will be in the subfolder "Plots"
if doCreateCirclesAnimation
    figuresToCreate='map'; % Can be also {'map','eigenVectors'}
    showFiguresWhileRuning='off'; % Whether to show or not the figures on screen

    for i=1:length(covStack)
        makecircle(covStack(i), doNormalize, dataFolder,figuresToCreate,showFiguresWhileRuning);
    end

    avconvCommand=['avconv -r 14 -i ' dataFolder '/Circle_plots/*.png -vcodec libx264 ' dataFolder filesep dataSetDescription '_' sprintf('%03d',round(covThreshold*100)) '.mp4'];
    avconvCommand=['avconv -r 14 -i ' dataFolder '/Covariance_animation_frames/%05d.png -vcodec libx264 -vf crop=860:870:40:0 ' dataFolder filesep dataSetDescription '_EOFs.mp4'];
    disp('To produce animation execute:')
    disp([' >' avconvCommand]);
end
%% CREATING COVARIANCE ANIMATIONS
if doCreateCovAnimation
    animateCovariances(covStack,covThreshold,dataFolder)

    avconvCommand=['avconv -r 14 -i ' dataFolder '/Covariance_animation_frames/%05d.png -vcodec libx264 -vf crop=860:870:40:0 ' dataFolder filesep dataSetDescription '_' sprintf('%03d',round(covThreshold*100)) '.mp4'];
    disp('To produce animation execute:')
    disp([' >' avconvCommand]);
end