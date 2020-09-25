function emaxevent()

% 2-compartiment model with stat events. Plasma concentration is kept
% within a given range by infusion rate: if drug concentration achieves
% macimum then decrease infusion; if minimum - infusion in increased.


% PK parameters
p.Rate = 1.0; % infusion rate
p.V = 6; % volume of distribution
p.k = 0.1; % elimination rate
p.Emax = 20; % max effect (hill equation)
p.E50 = 10; % concentration linked to 50% max effect

options = odeset('Events', @fevents);

% initial values
tspan = [0 150]; % max time domain
c0 = 0; % initial value for the dependent variable
ie = 99; % avid empty events vector

while ~isempty(ie)
    % Call ode solver
    [t, c, te, ye, ie] = ode45(@derivatives, tspan, c0, options, p);
    % generate plot
    plot(t, c(:,1), '-', 'Color', 'k', 'LineWidth', 3)
    hold on;
    plot(t, hillEffect(c(:,1), 0, p.Emax, p.E50, 1), 'Color', 'b', 'LineWidth', 3)
    
    tspan = [te(end) 500];
    c0 = ye(end,1);
    R1 = p.Rate;
    % check if an event is detected
    switch ie(end)
        case 1 % drug effect E has increased to 2 - stop infusion
            p.Rate = 0;
        case 2 % drug effect E has decreated to 1 - restart infusion
            p.Rate = 1;
        case 3
            ie = []; % end of simulation - no more events
    end

    R2 = p.Rate;
    plot([t(end) t(end)], [R1 R2], 'Color', 'r', 'LineWidth',2)
    
end

grid minor;
grid on;
title('')
xlabel('Time [h]');
ylabel({'Drug effect', 'Concentration [mg/L]', 'Infusion Rate [mg/h]'});
legend('rate', 'conc', 'effect')

%save plot
print('-r900', '-dtiff', 'emaxevent')
end

function dcdt = derivatives(~,c,p)
% model

dcdt = p.Rate/p.V - p.k*c(1);
end

function [value, isterminal, direction] = fevents(t,y,p)
% fevents specify events for ode model.

% event 1
% % drug effect has increased to the upper limit
value(1,1) = hillEffect(y(1), 0, p.Emax, p.E50, 1) - 0.5;
isterminal(1,1) = 1; % stop the integrator
direction(1,1) = 1; % positive direction only

% event 2
% drug effect E has decreased to the lower limit
value(2,1) = hillEffect(y(1), 0, p.Emax, p.E50, 1) - 0.3; % end of observation
isterminal(2,1) = 1; % stop the integrator
direction(2,1) = -1; % positiive direction only

% event 3
% end of time integration
value(3,1) = t - 25; % detect concentration = max while growing
isterminal(3,1) = 1; % stop the integrator
direction(3,1) = 1; % positiive direction only
end

function E = hillEffect(c, E0, Emax, EC50, n)
    E = E0 + Emax.*c.^n./(EC50.^n+c.^n);
end