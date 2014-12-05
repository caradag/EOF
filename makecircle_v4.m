function results = makecircle_v4(cov_data, yr, wk, norm, workpath, plots, showFigures)
%% v4 has a revised normalizing factor for the amplitude series and only 
%  plots the pressure series belonging to the biggest component of the 
%  eigenvector.
% Plotting the e/vectors with eigenvalues of the covariance matrix that 
% satisfy the lambda test, then plotting the components of the eigenvector 
% above thres as circles (red for positive and blue otherwise; circle size 
% represents magnitude) on the relative position of sensors map. 
% Inputs: 
%    cov_data - string, name of the m.file containg the covariance 
%           matrix and their corresponding sensor's ID.  
%           Eg: '5days_int_unnorm_cov_data_2008_1.mat'
%    yr -   scalar, the year of the cov_data.
%    wk -   scalar, the week of the cov_data.
%    norm - Boolean, user-specified. True for normalized data. 
%           Default is false.
%    workpath - path were figures will be saved to and were cov_data is stored
%    plots - Cell vector of strings with the lists of plot names to be
%           produced.
%           Default is all, i.e. plots={'eigenVectors','map'};
%    showFigures - 'on' or 'off' set the visibility of the figures.
%           this is if they will be displayed and saved or just saved to disk
%           Default is 'on'
%
% Output: 
%    results - string, irrelevant. See saved plots.
%

% Setting user-specified inputs to default if not specified
if nargin < 4 || isempty(norm);
    norm=false;
end

if nargin < 5 || isempty(workpath);
    workpath=[pwd filesep];
end

if nargin < 6 || isempty(plots);
    plots={'eigenVectors','map'};
end
if nargin < 7 || ~any(strcmp(showFigures,{'on','off'}));
    showFigures='on';
end

persistent data sensors northing easting

if isempty(data)
    tic;
    fprintf('Loading data file...\n');
    data = load('data 2014 v5 good only.mat');
    fprintf('Loaded in %.1f seconds.\n',toc);
    sensors = fieldnames(data);
end

% Loading data and setting up some useful variables
% Loaded variables are (as in Dec 4th 2014):
%   gridpos           192x1              14984  cell   (unused)             
%   positions         192x1             163584  cell                
%   sensor_count        1x1                  8  double              
%   sensors           192x1              13824  cell 
if isempty(northing) || isempty(easting)
    tic;
    fprintf('Loading sensor locations file...\n');    
    load('location of sensors.mat');
    fprintf('Loaded in %.1f seconds.\n',toc);
    sensor_count = length(fieldnames(data));
    northing = nan(sensor_count,1); easting = nan(sensor_count,1);
    for ii=1:sensor_count;
        northing(ii)=positions{ii}.north;
        easting(ii)=positions{ii}.east;
        % Changing to nominal coordinates if coordinates are not available
        if isnan(northing(ii))||isnan(easting(ii));
            northing(ii)=positions{ii}.nominal_north;
            easting(ii)=positions{ii}.nominal_east;
        end
    end    
end

load([workpath cov_data]);
[V,D] = eig(cov_clean); % D is a diagonal matrix of e/values and V is a matrix  
                        % whose columns are the corresponding eigenvectors
D2 = diag(D); % D2 is a column vector of e/values

% Getting the e/values (eigval2) and e/vectors (V2) that satisfies the lambda test
lambdabar = sum(D2)/length(sensors_clean);
ind = D2>lambdabar;
eigval = D2(ind); V2 = V(:,ind); % eigval is a column vector containing the 
                                 % eigenvalues

% Normalizing the e/vector (V2)
[m,n] = size(V2);
eigvec_mat = V2./repmat(max(abs(V2),[],1),m,1); % eigvec_mat is a matrix whose 
                                                % columns are the eigenvectors

% Setting up some useful variables for the ampltiude time series
dt = 2/60/24;
% A_mat is a matrix whose columns are going to be the amplitude series
m2 = length(datenum(start_time):dt:datenum(end_time));
A_mat = zeros(m2,n); 

startTime=datenum(start_time);
endTime=datenum(end_time);

% Plotting
% Looping through the eigenvectors
for jj=1:n;
    % Setting/resetting A to a column vector of zeros
    A = zeros(m2,1);
    
    % Changing the sign of the component values of the e/vectors if the most extreme
    % component value is negative
    if max(eigvec_mat(:,jj))<1;
        eigvec_mat(:,jj) = -eigvec_mat(:,jj);
    end
    eigvec = eigvec_mat(:,jj); % eigvec is the eigenvector jj
    
    if any(strcmp(plots,'eigenVectors'))
        figHandle=figure('Visible', showFigures);
        plot(eigvec,'bx--');
        xlabel('Component'); ylabel('Magnitude'); 
        filename = [num2str(yr) '_' num2str(wk) '_Eigvec ' num2str(jj)];
        title(['Eigenvector ' num2str(jj) ' (' num2str(yr) '-' num2str(wk) ' ' datestr(startTime,'mmm-dd')  ' to ' datestr(endTime,'mmm-dd') ')'],'Color','w');
        saveas(gcf, [workpath 'Plots/' filename], 'png');
        delete(figHandle);
    end 
    
    if any(strcmp(plots,'map'))
        % Creating a red/blue transparent circle around the sensors, plotting the 
        % pressure series belonging to the biggest eigenvector component and the
        % amplitude series

        % Looping through the sensors
        figHandle=figure('Name',['Map #' num2str(jj)],'Visible', showFigures);
        subplot(3,5,[1:3 6:8 11:13]);
        box on
        set(gca,'Color',[.99 .99 .99],'xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99],'FontSize',9,'XTickLabel',{},'YTickLabel',{});
        hold on
        % Plotting the position of the sensors
        plot((easting-min(easting)),(northing-min(northing)),'kx','MarkerSize',4);
        title(['Map ' num2str(jj) ' (' num2str(yr) '-' num2str(wk) ', ' datestr(startTime,'mmm-dd')  ' to ' datestr(endTime,'mmm-dd') ')'],'Color',[.99 .99 .99],'FontSize',15);
        xlim([0 max(easting)-min(easting)]);
        ylim([0 max(northing)-min(northing)]);
        
    end
    for kk=1:m;
        % Assigning data as inputs for extract_v4
        t_in = data.(sensors_clean{kk}).time.serialtime;
        p_in = data.(sensors_clean{kk}).pressure;
        p_in = p_in{1};
        tlim = [min(t_in) max(t_in)];
        if datenum(start_time)<tlim(1) || datenum(end_time)>tlim(2);
            % Making all the entries in p_out zero if datenum(start_time)<tlim(1) or 
            % datenum(end_time)>tlim(2)
            t_out = (datenum(start_time):dt:datenum(end_time)).';
            p_out = zeros(length(t_out),1);
        else
            % Calling extract_v6
            [p_out, t_out] = extract_v6(datenum(start_time), datenum(end_time),t_in, p_in, norm, dt);
            if any(isnan(p_out));
                p_out = zeros(length(t_out),1);
            end
        end
        c = sqrt(sum(eigvec.^2)); % Normalizing factor
        A = A + p_out.*eigvec(kk)./c;

       if any(strcmp(plots,'map'))
            % Plotting the pressure series belonging to the biggest eigenvector component
            if eigvec(kk)==max(eigvec);
                % Setting up variables for patch
                x1 = datenum(start_time);
                x2 = datenum(end_time);
                y1 = min(p_in);
                y2 = max(p_in);

                % Back to plotting                
                subplot(3,5,[4 5]);
                set(gca,'Color',[.99 .99 .99],'xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99],'FontSize',10);
                % Plotting the the whole time serie decimated to 300 samples
                hold on
                box on
                patch([x1 x2 x2 x1],[y1 y1 y2 y2]/(1000*9.8),'yellow');
                decimationFac=round(length(t_in)/300);
                decimatedIndexes=1:decimationFac:length(t_in);
                plot(t_in(decimatedIndexes), p_in(decimatedIndexes)/(1000*9.8),'b.','MarkerSize',3);
                title(['Press. ' sensors_clean{kk}],'Color',[.99 .99 .99],'FontSize',15);
                datetick('x','yyyy');
                
                subplot(3,5,[9 10]);
                set(gca,'Color',[.99 .99 .99],'xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99],'FontSize',10);
                
                % Plotting the data inside the window at full resolution
                hold on
                box on
                inWindow= t_in>=datenum(start_time) & t_in<=datenum(end_time);
                plot(t_in(inWindow), p_in(inWindow)/(1000*9.8),'b.','MarkerSize',3);
                title('Window press.','Color',[.99 .99 .99],'FontSize',10);
                datetick('x','mm-dd');
                
            end

            % Creating red/blue transparent circles
            subplot(3,5,[1:3 6:8 11:13]);
            hold on
            ind2 = strcmp(sensors_clean{kk},sensors);
            if eigvec(kk)>=0;
                c = 'r';
            else
                c = 'b';
            end
            xc = easting(ind2)-min(easting);
            yc = northing(ind2)-min(northing);
            r = eigvec(kk)*30;
            x = r*sin(-pi:0.1*pi:pi) + xc;
            y = r*cos(-pi:0.1*pi:pi) + yc;
            h=fill(x, y, c, 'FaceAlpha', 0.4);
            uistack(h,'top');
            text(xc, yc, num2str(kk));
       end
    end

    % Setting the column of A_mat to A
    A_mat(:,jj) = A;

    if any(strcmp(plots,'map'))    
%         subplot(3,5,[1:3 6:8 11:13]);
%         set(gca,'Position',get(gca,'OuterPosition'));
        
        % Plotting the amplitude against time graph (subplot(2,1,2))
        subplot(3,5,[14:15]);
        set(gca,'xcolor',[.99 .99 .99],'ycolor',[.99 .99 .99],'Color',[.99 .99 .99],'FontSize',9)
        
        hold on
        box on
        plot(t_out, A_mat(:,jj)/(9.8*1000),'b.','MarkerSize',3);
        xlabel('Time/days');
        title('Amplitude vs time graph','Color',[.99 .99 .99],'FontSize',9);
        datetick('x','mm-dd');
        
        % Saving the plot
        filename = [num2str(yr) '_' num2str(wk) '_Map ' num2str(jj)];
        disp(filename);
        
        saveas(gcf, [workpath 'Plots/' filename], 'png');
        %print(gcf,'-dpng','-r600',[workpath 'Plots/' filename '.png'])
        if strcmp(showFigures,'on')
            input('Continue?')
            pause(3);
        end
        delete(figHandle);
    end
end


% Saving the variables
filename = ['E_val, e_vec, amp and time' ' (' cov_data(1:end-4) ')'];
save([workpath 'Eig and amp/' filename], 'eigval', 'eigvec_mat', 'A_mat','start_time', 'end_time');

close ALL
results = 'See saved plots';
end