function results = makecircle_v4(cov_data, yr, wk, norm)
%% v4 has a revised normalizing factor for the amplitude series and only 
%  plots the pressure series belonging to the biggest component of the 
%  eigenvector.
% Plotting the e/vectors with eigenvalues of the covariance matrix that 
% satisfy the lambda test, then plotting the components of the eigenvector 
% above thres as circles (red for positive and blue otherwise; circle size 
% represents magnitude) on the relative position of sensors map. 
% Inputs: cov_data - string, name of the m.file containg the covariance 
%                    matrix and their corresponding sensor's ID.  
%                    Eg: '5days_int_unnorm_cov_data_2008_1.mat'
%         yr - scalar, the year of the cov_data.
%         wk - scalar, the week of the cov_data.
%         norm - Boolean, user-specified. True for normalized data. 
%                Default is false.
% Output: results - string, irrelevant. See saved plots.

% Setting user-specified inputs to default if not specified
if nargin < 4 || isempty(norm);
    norm=false;
end

data = load('data 2014 v5 good only.mat');
sensors = fieldnames(data);
sensor_count = length(fieldnames(data));

% Loading data and setting up some useful variables
load('location of sensors.mat');
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

load(cov_data);
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
    plot(eigvec,'bx--');
    xlabel('Component'); ylabel('Magnitude'); 
    filename = [num2str(yr) '_' num2str(wk) '_Eigvec ' num2str(jj)];
    title(['Eigenvector ' num2str(jj) ' (' start_time ' to ' end_time ')']);
    saveas(gcf, ['Plots/' filename], 'fig');
    clf(gcf);
    
    % Creating a red/blue transparent circle around the sensors, plotting the 
    % pressure series belonging to the biggest eigenvector component and the
    % amplitude series
    
    % Looping through the sensors
    figure(1);
    for kk=1:m;
        % Assigning data as inputs for extract_v4
        t_in = eval(['data.' sensors_clean{kk} '.time.serialtime']);
        p_in = eval(['data.' sensors_clean{kk} '.pressure']);
        p_in = p_in{1};
        tlim = [min(t_in) max(t_in)];
        if datenum(start_time)<tlim(1) || datenum(end_time)>tlim(2);
            % Making all the entries in p_out zero if datenum(start_time)<tlim(1) or 
            % datenum(end_time)>tlim(2)
            t_out = (datenum(start_time):dt:datenum(end_time)).';
            p_out = zeros(length(t_out),1);
        else
            % Calling extract_v6
            [p_out, t_out] = extract_v6(datenum(start_time), datenum(end_time),...
                t_in, p_in, norm, dt);
            if any(isnan(p_out));
                p_out = zeros(length(t_out),1);
            end
        end
        c = sqrt(sum(eigvec.^2)); % Normalizing factor
        A = A + p_out.*eigvec(kk)./c;
        
        % Plotting the pressure series belonging to the biggest eigenvector component
        if eigvec(kk)==max(eigvec);
            % Setting up variables for patch
            x1 = datenum(start_time)-datenum('1-Jan-2008');
            x2 = datenum(end_time)-datenum('1-Jan-2008');
            y1 = min(p_in);
            y2 = max(p_in);
            
            % Back to plotting
            subplot(2,5,[4:5]);
            patch([x1 x2 x2 x1],[y1 y1 y2 y2],'yellow');
            hold on
            plot((t_in-datenum('1-Jan-2008')), p_in,'b*');
            xlabel('Time/days');
            ylabel('Pressure')
            title(['Pres. Series, sensor ' sensors_clean{kk} ' (' num2str(kk) ')']);
            hold off
        end
        
        % Creating red/blue transparent circles
        subplot(2,5,[1:3 6:8]);
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
        fill(x, y, c, 'FaceAlpha', 0.4);
        text(xc, yc, num2str(kk));
        hold on
    end
    
    % Setting the column of A_mat to A
    A_mat(:,jj) = A;
    
    % Plotting the position of the sensors
    plot((easting-min(easting)),(northing-min(northing)),'kx');
    xlabel('Easting/m');
    ylabel('Northing/m');
    title(['Map ' num2str(jj) ' (' start_time ' to ' end_time ')']);
    hold off

    % Plotting the amplitude against time graph (subplot(2,1,2))
    subplot(2,5,[9:10]);
    plot(t_out-min(t_out), A_mat(:,jj));
    xlabel('Time/days');
    ylabel('Amplitude');
    title('Amplitude vs time graph');
    
    % Saving the plot
    filename = [num2str(yr) '_' num2str(wk) '_Map ' num2str(jj)];
    saveas(gcf, ['Plots/' filename], 'fig');
    clf(gcf);
end

% Saving the variables
filename = ['E_val, e_vec, amp and time' ' (' cov_data(1:end-4) ')'];
save(['Eig and amp/' filename], 'eigval', 'eigvec_mat', 'A_mat',...
    'start_time', 'end_time');

close ALL
results = 'See saved plots';
end