function infus2com()

% 2-compartiment model with infusion. In the treatment, doses are given as
% infusions to each subject. A graph with time courses of plasma
% concentration and infusion is generated. The function has no arguments.

p.Rate = 50; % initial infusion rate = 100/2 [mg/h]
p.CL = 1.38e01; % central clearance.
p.V1 = 1.48e01; % volume of distribution in central compartment.
p.Q = 2.17e00; % inter-compartmental clearance.
p.V2 = 4.23e00; % volume of distribution peripheral compartment.
p.k = p.CL/p.V1; % rate constant of elimination
p.k12 = p.Q/p.V1; % rate constant from central to peripheral.
p.k21 = p.Q/p.V2; % rate constant from peripheral to central

options = odeset('Events', @fevents, 'RelTol', 1e-6, 'AbsTol',...
    [1e-6 1e-6]);

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
        case 1 % event: stop infusion
            p.Rate = 0;
        case 2 % event: end of observation
            ie = [];
        otherwise
            % no action needed
    end
end

% plot
title('2-compartment model with infusion')
xlabel('Time [h]');
ylabel('Concentration [mg/L]');
text('units','inch','position',[1.9 2.7],...
    'fontsize',11, 'interpreter', 'latex', 'string',...
    ['$${dc_1 \over dt} = {r\over V_1} - (k + k_{12})\cdot c_1 + '...
    'k_{21}\cdot {V_2 \over V_1} \cdot c_2 $$'])
text('units', 'inch', 'position',[1.9 2.2],...
    'fontsize',11, 'interpreter','latex','string',...
    ['$${dc_2 \over dt} = k_{12} \cdot {V_1 \over V_2}'...
    '\cdot c_1 - k_{21} \cdot c_2 \ $$'])
legend('central compartment', 'peripheral compartment')

%save plot
print('-r900', '-dtiff', 'infus2com')
end

function dcdt = derivatives(~,c,p)
% model

dcdt = [p.Rate/p.V1 - (p.k + p.k12)*c(1) + p.k21*p.V2/p.V1*c(2)
    p.k12*p.V1/p.V2*c(1) - p.k21*c(2)];
end

function [value, isterminal, direction] = fevents(t,~,~)
% fevents specify events for ode model.

% 1. Event: Locate the time to stop infusion
value(1,1) = t - 2; % detect t=2 (stop time)
isterminal(1,1) = 1; % stop the integrator
direction(1,1) = 1; % positive direction only

% 2. Event: Locate the time to stop the simulation
value(2,1) = t - 15; % detect t=15 (stop time)
isterminal(2,1) = 1; % stop the integrator
direction(2,1) = 1; % positiive direction only
end