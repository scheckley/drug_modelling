function PK2IV1
% solev 2-compartment model with bolus input and 1st order elimination. The
% mass option is applied.

clear all;
close all;
clc;

p.V = 5.0; % L
p.k10 = 1.0; % mg/h
p.k12 = 0.8; % mg/h
p.k21 = 0.3; % mg/h

dose = 1000; % mg

timeEnd = 10; % h

options = odeset;
options.Mass = [1 0 0; 0 0 0; 0 0 1];

[t, y] = ode15s(@derivalgeb, [0, timeEnd], [dose; dose/p.V;0],...
    options, p);

plot(t, y(:,1), '-r', t, y(:,2),'b-', t,y(:,3), 'g-', 'lineWidth',2)
legend('Central Compartment: drug amount [mg]',...
    'Central Compartment: Drug concentrations [mg/L]',...
    'Peripheral Compartment: Drug amounts [mg]')
xlabel('Time (h)')
ylabel('Value')
figure()
plot(t,y(:,2), 'b-', 'lineWidth', 2)
legend('Drug concentrations [mg/L]')
xlabel('Time [h]')
ylabel('Value')
%print('-dtiff','-r900', 'PK2IV1.tif')

end

function dydt = derivalgeb(~,y,p)
a1 = y(1);
c1 = y(2);
a2 = y(3);
da1 = -p.k10 * a1 - p.k12 * a1 + p.k21 * a2;
dc1 = c1 - a1/p.V;
da2 = p.k12 * a1 - p.k21 * a2;
dydt = [da1; dc1; da2];

end



