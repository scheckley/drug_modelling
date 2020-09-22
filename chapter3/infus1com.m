function infus1com

% one compartment model with infusion
% infus1com computes plasma concentration versus time course and creates
% plots for a one-compartment model. Drug doses are applied as infusions.
% The function has no parameters.

% PK parameters

p.Dose = 200; %dose [mg]
p.tStop = 2; % infusion time stop [h]
p.R = p.Dose/p.tStop; % infusion rate [mg/h]
p.V = 12.6; % volume of distribution [L]
p.Cl = 2.16; % clearance [L/h]


options = odeset('Events', @fevents, 'RelTol', 1e-3, 'AbsTol', 1e-3);

ie = 99; % any non-empty value to avoid 't=0' as event
tspan = [0,100]; % integration time interval
c0 = 0; % initial concentration value

% prepare data for plotting

timePoints = [];
concPoints = [];
ratePoints = [];

% Call ode solver and check events

while ~isempty(ie)
    [t, c, te, ce, ie] = ode45(@derivatives, tspan, c0, options, p);
    % collect data for plotting
    timePoints = [timePoints; t];
    concPoints = [concPoints; c];
    ratePoints = [ratePoints; p.R*ones(length(t),1)];
    tspan = [te(end) 100];
    c0 = ce(end,1);
    %check if an event is detected
    switch ie(end)
        case 1 % event: stop infusion
            p.R = 0;
        case 2 % event: end of observation
            ie = [];
        otherwise
            % no action needed - placeholder for other events
    end
end

% generate plots
subplot(3,1,1)
plot(timePoints, ratePoints, '-k', 'LineWidth', 2);
hold on
title({' ', '\fontsize{13}One-Compartment Model with Infusion', ''})
ylim([0 110])
ylabel('\fontsize{13}Rate [mg/h]')
subplot(3,1,2:3);
plot(timePoints, concPoints, '-k', 'LineWidth',2)
hold on;
grid minor;
grid on;
xlabel('\fontsize{13}Time [hour]')
ylabel('\fontsize{13}Concentration [mg/L]')

% save plots
print('-r900', '-dtiff', 'infus1com')
end

function dcdt = derivatives(t, c, p)
% the model
dcdt = r(t, p)/p.V - p.Cl/p.V*c;
end

function rt = r(~, p)

% r specifies the infusion rate
% rt = r(t,p) returns infusion rate 'rt' at time 't', and with parameters
% 'p'.

% infusion (depends on event)
rt = p.R;
end

function [value, isTerminal, direction] = fevents(t, ~, p)
%[value, isterminal, direction] = fevents(t,y,p) detects the exact time
%points of events specified as a function of time 't', vector of
%independent variables 'y', and a structure of some (specific) other
%parameters 'p'.

% Output parameters: 
% 'value' - value investigated as event
% 'isterminal' = 0/1, where 0 stands for 'don't stop', 1 for 'stop
% integration'.
% 'direction' = -1/0/1 where -1 means 'negative direction
% only', 0 - 'both direction', 1 - 'positive direction only'

% 1. Event: Locate the time to stop infusion
value(1,1) = t - p.tStop; % Detect t=2 (infusion stop time)
isTerminal(1,1) = 1; % Interrupt integration, stop infusion
direction(1,1) = 1; % Positive direction only

% 2. Event: Locate the time to stop observation/integration
value(2,1) = t - 15; % Detect t=15 (end of observation)
isTerminal(2,1) = 1; % Stop integration
direction(2,1) = 1; % Positive direction only
end
