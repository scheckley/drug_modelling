function infus2comstev()

% 2-compartiment model with stat events. Plasma concentration is kept
% within a given range by infusion rate: if drug concentration achieves
% macimum then decrease infusion; if minimum - infusion in increased.


% PK parameters
CL = 1.38e01; % central clearance.
V1 = 1.48e01; % volume of distribution in central compartment.
Q = 2.17e00; % inter-compartmental clearance.
V2 = 4.23e00; % volume of distribution peripheral compartment.
Rate = 50; % infusion rate

% p structure with PK parameters available in several functions
p.Rate = Rate;
p.V1 = V1;
p.V2 = V2;
p.k = CL/V1;
p.k12 = Q/V1;
p.k21 = Q/V2;

options = odeset('Events', @fevents);

% initial values
tspan = [0 50]; % max time domain
c0 = [0 0]; % initial value for the dependent variable
ie = 99; % avid empty events vector

while ~isempty(ie)
    % Call ode solver
    [t, c, te, ye, ie] = ode45(@derivatives, tspan, c0, options, p);
    % generate plot
    plot(t, c(:,1), '-', 'Color', 'r', 'LineWidth', 2)
    hold on;
    plot(t, c(:,2), '-.', 'Color', 'b', 'LineWidth', 2)
    
    tspan = [te(end) 50];
    c0 = [ye(end,1) ye(end,2)];
    % check if an event is detected
    switch ie(end)
        case 1 % event: stop infusion after 2 hours
            p.Rate = 0;
        case 2 % event: stop process definitely
            ie = [];
        case 3 % reduce infusion if conc. has grown to allowed max
            p.Rate = Rate/10;
        case 4 % allow full infusion rate if conc. has fallen to min
            p.Rate = Rate;
        otherwise
            messageID = 'infus2comstev:UnknownEvent';
            messageStr='Unknown event ''%s''.';
            error(messageID, messageStr, num2str(ie(end)));
    end
end

% plot
hold off;
title('2-compartment model with controlled infusion rate')
xlabel('Time [h]');
ylabel('Concentration [mg/L]');
legend('central compartment', 'peripheral compartment')

%save plot
print('-r900', '-dtiff', 'infus2com')
end

function dcdt = derivatives(~,c,p)
% model

dcdt = [p.Rate/p.V1 - (p.k + p.k12)*c(1) + p.k21*p.V2/p.V1*c(2)
    p.k12*p.V1/p.V2*c(1) - p.k21*c(2)];
end

function [value, isterminal, direction] = fevents(t,y,~)
% fevents specify events for ode model.

value(1,1) = t - 2; % max time for first infusion
isterminal(1,1) = 1; % stop the integrator
direction(1,1) = 1; % positive direction only

value(2,1) = t - 30; % end of observation
isterminal(2,1) = 1; % stop the integrator
direction(2,1) = 1; % positiive direction only

value(3,1) = y(1) - 2; % detect concentration = max while growing
isterminal(3,1) = 1; % stop the integrator
direction(3,1) = 1; % positiive direction only

value(4,1) = y(1) - 0.5; % detect concentration = min while falling
isterminal(4,1) = 1; % stop the integrator
direction(4,1) = -1; % negative direction only
end