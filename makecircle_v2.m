function results = makecircle_v2(cov_data, yr, wk, thres1, thres2, norm)
%% v2 is used in cleanup; new plots have "(2)" at the end of their names. 
% Plotting the e/vectors with eigenvalues of the covariance matrix above thres1, 
% then plotting the components of the eigenvector above thres2 as circles (red for  
% positive and blue otherwise; circle size represents magnitude) on the relative 
% position of sensors map. 
% Inputs: cov_data - string, name of the m.file containg the covariance matrix and 
%                    their corresponding sensor's ID. 
%                    Eg: '5days_int_unnorm_cov_data_2008_1.mat'
%         yr - scalar, the year of the cov_data.
%         wk - scalar, the week of the cov_data.
%         thres1 - scalar (0 to 1.0), user-specified. The percentage of cumulative 
%                  energy (sum of e/values) that decides the e/vectors that are   
%                  to be included in the plots. Default is 0.9.  
%         thres2 - scalar (0 to 1.0), user-specified. The miniumum absolute value   
%                  of the eigenvector component for which the pressure series will be 
%                  plotted. Default is 0.5.  
%         norm - Boolean, user-specified. True for normalized data. Default is false.
% Output: results - string, irrelevant. See saved plots.

% Setting user-specified inputs to default if not specified
if nargin < 4 || isempty(thres1);
    thres1=0.9; 
end

if nargin < 5 || isempty(thres2);
    thres2=0.5; 
end

if nargin < 6 || isempty(norm);
    norm=false;
end

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

load('clean data v6.mat');

load(cov_data);
[V,D] = eig(cov_clean); % D is a diagonal matrix of e/values and V is a matrix whose 
                        % columns are the corresponding eigenvectors
D2 = diag(D); % D2 is a column vector of e/values
totalener = sum(D2);

% Getting the e/values (eigval) and e/vectors (V2) that satisfies thres1
ener = 0; count = 0;
while ener<thres1*totalener;
    ener = ener+D2(end-count);
    count = count+1;
end
eigval = D2(end-count+1:end); 
V2 = V(:,length(sensors_clean)-count+1:length(sensors_clean));

% Normalizing the e/vector (V2)
[m,n] = size(V2);
eigvec_mat = V2./repmat(max(abs(V2),[],1),m,1); % eigvec_mat is a matrix whose 
                                                % columns are the eigenvectors

% Setting up some useful variables for the ampltiude time series
dt = 2/60/24;
% A_mat is a matrix whose columns are going to be the amplitude series
m2 = length(datenum(start_time):dt:datenum(end_time));
A_mat=zeros(m2,n); 

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
    filename = [num2str(yr) '_' num2str(wk) '_Eigvec ' num2str(jj) ' (2)'];
    title(['Eigenvector ' num2str(jj) ' (' start_time ' to ' end_time ')']);
    saveas(gcf, filename, 'fig');
    clf(gcf);
    
    % Creating a red/blue transparent circle around the sensors (figure(1), 
    % subplot(2,1,1)) and plotting the pressure series if val(kk)>thres2
    
    % Looping through the sensors
    figure(1);
    subplot(2,1,1);
    for kk=1:m;
        % Assigning data as inputs for extract_v4
        t_in = eval(['data.' sensors_clean{kk} '.time.serialtime']);
        p_in = eval(['data.' sensors_clean{kk} '.pressure']);
        p_in = p_in{1};
        tlim = eval(['data.' sensors_clean{kk} '.limits.time']);
        if datenum(start_time)<tlim(1) || datenum(end_time)>tlim(2);
            % Making all the entries in p_out zero if datenum(start_time)<tlim(1) or 
            % datenum(end_time)>tlim(2)
            t_out = (datenum(start_time):dt:datenum(end_time)).';
            p_out = zeros(length(t_out),1);
        else
            % Calling extract_v4
            [p_out, t_out] = extract_v3(datenum(start_time), datenum(end_time),...
                t_in, p_in, dt, norm);
        end
        c = sum(eigvec); % Normalizing factor
        A = A + p_out.*eigvec(kk)./c;
        
        % Plotting and saving the pressure series if thres2 is exceeded
        if abs(eigvec(kk))>thres2;
            figure(2);
            plot(t_out-min(t_out), p_out);
            xlabel('Time/days');
            ylabel('Pressure')
            filename = ['Pressure Series, sensor ' num2str(kk) ' (' start_time...
                ' to ' end_time ', Map ' num2str(jj) ') (2)'];
            title(filename);
            saveas(gcf, filename, 'fig');
            clf(gcf); figure(1);
        end
        
        % Creating red/blue transparent circles
        ind2 = strcmp(sensors_clean{kk},sensors);
        if eigvec(kk)>=0;
            c = 'r';
        else
            c = 'b';
        end
        xc = easting(ind2);
        yc = northing(ind2);
        r = eigvec(kk)*10;
        x = r*sin(-pi:0.1*pi:pi) + xc;
        y = r*cos(-pi:0.1*pi:pi) + yc;
        fill(x, y, c, 'FaceAlpha', 0.4);
        text(xc, yc, num2str(kk));
        hold on
    end
    
    % Setting the column of A_mat to A
    A_mat(:,jj) = A;
    
    % Plotting the position of the sensors
    plot(easting,northing,'kx');
    xlabel('Easting/m');
    ylabel('Northing/m');
    title(['Map ' num2str(jj) ' (' start_time ' to ' end_time ')']);
    hold off

    % Plotting the amplitude against time graph (subplot(2,1,2))
    subplot(2,1,2);
    plot(t_out-min(t_out), A);
    xlabel('Time/days');
    ylabel('Amplitude');
    title('Amplitude vs time graph');
    
    % Saving the plot
    filename = [num2str(yr) '_' num2str(wk) '_Map ' num2str(jj) ' (2)'];
    saveas(gcf, filename, 'fig');
    clf(gcf);
end
% Saving the variables
filename = ['E_val, e_vec, amp and time' ' (' cov_data(1:end-4) ') (2)'];
save(filename, 'eigval', 'eigvec_mat', 'A_mat', 'start_time', 'end_time');

close ALL
results = 'See saved plots';
end