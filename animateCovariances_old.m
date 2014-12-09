%Script that can animate covariance data as in the old kevin's folder structure

% clear all;
% covStack=[];

yearstart = 2008; 
yearend = 2014;

startWeek=1;
endWeek=73;

windowLength=5; % To process data of 5 days windows
windowOffsets=5; % Number of offset in the window for example:
                 %      5 => 5 days window with 1 day time steps
                 %      3 => 15 days window with 5 day time steps
normData = true; % true to process normalized data, false for unnormalized

if isempty(covStack)
    % Looping through the years
    for year=yearstart:yearend
        % Looping through the 'weeks'
        for week=startWeek:endWeek
            %looping trough folders identified by Day Of Window number
            for dow=1:windowOffsets
                if normData
                    folderName=sprintf('Results/Normalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
                    cov_data = sprintf('%ddays_int_norm_cov_data_%d_%d.mat',windowLength,year,week);
                else
                    folderName=sprintf('Results/Unnormalized press, %d days int, 2008 to 2014 (%d)/',windowLength,dow);
                    cov_data = sprintf('%ddays_int_unnorm_cov_data_%d_%d.mat',windowLength,year,week);
                end
                if exist([folderName cov_data],'file')
                    load([folderName cov_data]);
                    n=length(sensors_clean);
                    if ~isempty(cov_clean) && n>1
                        covStack(end+1).cov=cov_clean(triu(true(n),1));
                        covStack(end).timeLims=[datenum(start_time) datenum(end_time)];                    
                        covStack(end).sensors=sensors_clean;                    
                        covStack(end).n=n;
                        covStack(end).year=year;
                        covStack(end).window=week;
                        covStack(end).dayOfWinwdow=dow;
                        fprintf('   Window %d out of %d DOW: %d (%d)\n',week,endWeek-startWeek+1,dow,year);
                    else
                        fprintf('   Window %d out of %d DOW: %d (%d) EMPTY\n',week,endWeek-startWeek+1,dow,year);
                    end
                end
            end
        end
    end
    
    toIncludePercent=98;
    covVals=sort(cat(1,covStack.cov));
    nVals=length(covVals);
    halfOut=round(nVals*(((100-toIncludePercent)/2)/100));
    minCov=covVals(halfOut);
    maxCov=covVals(end-halfOut);
    fprintf('%d%% of covariance samples between %.2g and %.2g\n',toIncludePercent,minCov,maxCov)

    metadata=load('data/data 2014 v5 good only_metadata.mat');
    sensors=fieldnames(metadata);
    minE=Inf;
    maxE=-Inf;
    minN=Inf;
    maxN=-Inf;
    for s=1:length(sensors)
        minE=min(minE,metadata.(sensors{s}).pos(1));
        maxE=max(maxE,metadata.(sensors{s}).pos(1));
        minN=min(minN,metadata.(sensors{s}).pos(2));
        maxN=max(maxN,metadata.(sensors{s}).pos(2));
    end
end


imageFile='/home/camilo/5_UBC/Data visualization GUI/Reference images/map1024.tif';
mapW=560;
mapH=560;
% mapFigure=figure('Name','Map overview','NumberTitle','off','Position',[2 60 mapW mapH],'Color',[1 1 1],'Visible','off');

mapFigure=figure('Name','Map overview','NumberTitle','off','Color',[1 1 1],'Visible','off','MenuBar','none');
pos=get(mapFigure,'Position');
pos(3)=mapW;
pos(4)=mapH;
set(mapFigure,'Position',pos);

baseImage=imread(imageFile);
baseImage=cat(3,flipud(baseImage(:,:,1)),flipud(baseImage(:,:,2)),flipud(baseImage(:,:,3)));
%reading reoreferenciation data
[pathstr, name, ~] = fileparts(imageFile);
[imageH, imageW, ~]=size(baseImage);
tfw=load([pathstr '/' name '.tfw']);
maxN=tfw(6);
minN=maxN+tfw(4)*(imageH-1);
minE=tfw(5);
maxE=minE+tfw(1)*(imageW-1);


mapAxes=axes();
set(mapAxes,'DataAspectRatio',[1 1 1],'YDir','normal');
set(mapAxes,'XTick',[],'YTick',[],'Units','pixels','XLim',[minE maxE],'YLim',[minN maxN],'Color',[1 1 1]);
axesPos=get(mapAxes,'Position');
reductionFactor=min(axesPos(3)/imageW,axesPos(4)/imageH);
baseImage=imresize(baseImage, reductionFactor);
[imageH, imageW, ~]=size(baseImage);

dotSize=11;
[dotE dotN]=meshgrid(1:dotSize,1:dotSize);
dotCenter=((dotSize-1)/2)+1;
dot=sqrt((dotE-dotCenter).^2 + (dotN-dotCenter).^2);


rangeE=linspace(minE,maxE,imageW);
rangeN=linspace(minN,maxN,imageH);
[E N]=meshgrid(rangeE,rangeN);


nWindows=length(covStack);
pairCount=0;
dotRadius=10;%Squeare of the radius of the dot
covThreshold=0.99;
semiAxisMinor=0.002;
for i=1:nWindows
    messageLength=fprintf('%06.2f%% %d/%d',100*pairCount/nVals,i,nWindows);
    if ~ishandle(mapFigure)
        break;
    end
    cla;
    hold on
    frame=baseImage;
    n=covStack(i).n;
    [cols rows]=meshgrid(1:n,1:n);
    cols=cols(triu(true(n),1));
    rows=rows(triu(true(n),1));
    nPairs=(n*n-n)/2;
    pairCount=pairCount+nPairs;
    px=[];
    py=[];
    maxFrameCov=0;
    plottedPairs=0;
    for pair=1:nPairs
        cov=covStack(i).cov(pair);
        if cov>0
            cov=cov/maxCov;
            color=logical([1 0 0]);
        else
            cov=cov/minCov;
            color=logical([0 0 1]);
        end
        cov=min(cov,1);
        maxFrameCov=max(maxFrameCov,cov);
        
        sensorID1=covStack(i).sensors{cols(pair)};
        sensorID2=covStack(i).sensors{rows(pair)};
        p1=metadata.(sensorID1).pos;
        p2=metadata.(sensorID2).pos;
        px=[px;p1(1);p2(1)];
        py=[py;p1(2);p2(2)];
        if cov<covThreshold
            continue
        end
        cov=(cov-covThreshold)/(1-covThreshold);
        plottedPairs=plottedPairs+1;
        pp = sqrt(sum((p1-p2).^2));
        [~,p1Col]=min(abs(rangeE-p1(1)));
        [~,p1Row]=min(abs(rangeN-p1(2)));
        [~,p2Col]=min(abs(rangeE-p2(1)));
        [~,p2Row]=min(abs(rangeN-p2(2)));
        maxWidth=ceil((pp*semiAxisMinor)/((maxE-minE)/imageW));
        margin=max(maxWidth,dotRadius);
        startRow=max(min(p1Row,p2Row)-margin,1);
        endRow=min(max(p1Row,p2Row)+margin,imageH);
        startCol=max(min(p1Col,p2Col)-margin,1);
        endCol=min(max(p1Col,p2Col)+margin,imageW);

        boxE=E(startRow:endRow,startCol:endCol);
        boxN=N(startRow:endRow,startCol:endCol);
        d1=sqrt((boxE-p1(1)).^2 + (boxN-p1(2)).^2);
        d2=sqrt((boxE-p2(1)).^2 + (boxN-p2(2)).^2);
        d=d1+d2;
        affectedInBox=d<pp*(1+semiAxisMinor) | d1<dotRadius | d2<dotRadius;
        vals=d(affectedInBox)/pp;
        affected=false(size(E));
        affected(startRow:endRow,startCol:endCol)=affectedInBox;
        
        vals=(1./vals);
        vals=(vals-min(vals))/(max(vals)-min(vals));
        vals=vals.^15;
        affected=cat(3,affected & color(1),affected & color(2),affected & color(3));
        frame(affected)=frame(affected)+uint8(vals*cov*255)/2;
    end
    image(rangeE,rangeN,frame);
    plot(px,py,'ok','MarkerSize',6,'MarkerFaceColor','g');

    title(sprintf('#%d: %d-%d (%d), %s (max. cov. %.2f, plotted pairs %d/%d)',i,covStack(i).year,covStack(i).window,covStack(i).dayOfWinwdow, datestr(covStack(i).timeLims(1)),maxFrameCov,plottedPairs,nPairs));
    drawnow;
    if ~ishandle(mapFigure)
        break;
    end    
    saveas(mapAxes,sprintf('frames/%05d.png',i),'png');
    fprintf('%c',8*ones(messageLength,1));    
end
