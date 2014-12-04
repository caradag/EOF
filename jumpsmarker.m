function [ind_jump, t_jump, delta_p_jump] = jumpsmarker(sensor, tinitial, tfinal, thres, threshold_type);
%% Marks the jump(s) and saves the pressure series plot. 1 Jan 2008 is set as day 0.
% Inputs: sensor - string; the sensor's ID. Eg: 'SO8P09'.
%         thres - the factor of the larger adjacent gradient for which the interval
%                 is considered a jump
% Output: ind_jump - index in the sensor time series corresponding to time
%                  stamp at end of jump
%         t_jump - time stamp at end of jump
%         delta_p_jump - size of jump

load('clean data v6.mat');
t = eval(['data.' sensor '.time.serialtime']);
p = eval(['data.' sensor '.pressure']);
p = p{1};

%keep track of indices in original time and data vector
ind0 = 1:length(t);

% Deleting the entries of p and t which has NaN as its t entry
p(isnan(t)) = []; ind0(isnan(t)) = []; t(isnan(t)) = [];
t(isnan(p)) = []; ind0(isnan(p)) = []; p(isnan(p)) = [];
p(t<tinital) = []; ind0(t<tinital) = []; t(t<tinitial) = [];
p(t>tfinal) = []; ind0(t>tfinal) = []; t(t>tfinal) = [];

% Getting the gradient
grad = abs(diff(p)./diff(t));

switch threshold_type
    case 'absolute'
        ind = grad>thres;
    case 'relative'
        ind = grad>thres*mean(grad);
    case 'change'   %fractional change in gradient exceeds threshold
        ind1 = grad>thres.*[grad(1);grad(1:end-1)];
        ind2 = grad>thres.*[grad(2:end);grad(end)];
        ind = ind1|ind2;
    otherwise
        error('not a recognized threshold type')
end

%output
t_jump = t([false; ind]);
delta_p_jump = p([false; ind])-p([ind; false]);
ind_jump = ind0([false; ind]);

% Setting up variables to mark the jumps on the plot
t2 = t([ind; false]);
p2 = p([ind; false]);
t3 = t([false; ind]);
p3 = p([false; ind]);

% Plotting
plot(t-datenum('1-Jan-2008'),p,'bx',t2-datenum('1-Jan-2008'),p2,'rx',...
    t3-datenum('1-Jan-2008'),p3,'rx');
xlabel('Time');
ylabel('Pressure');
title(['Pressure Series, Sensor ' sensor]);
end