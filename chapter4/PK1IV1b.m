function PK1IV1b

% simultaneous calculation of drug amounts and concentrations using a
% one-component PK model with IV input and first-order elimination. The 2nd
% equation defiend in PK1IVb is treated an an algebreic equations using the
% Mass option in odeset.

clear all;
clc;
close all;

p.V = 5.0; % L
p.k10 = 1.0; % mg/h
dose = 100; % mg
timeEnd = 10; % h

options = odeset;
options.Mass = [1 0; 0 0];

[t, y] = ode15s(@derivalgeb, [0 timeEnd], [dose; dose/p.V], options, p);

a = y(:,1);
c = y(:,2);

plot(t, a, '-r', t, c, 'b-', 'lineWidth', 2)
legend('Drug amount [mg]', 'Drug concentrations [mg/L]')
xlabel('Time [h]')
ylabel('Value')

end

function dydt = derivalgeb(~, y, p)
% derivalgeb compute the right hand side of the DAW

a1 = y(1);
c1 = y(2);
da1 = -p.k10 * a1; % ode
dc1 = c1 - a1/p.V; % algebreic equation: dc1 = 0
dydt = [da1; dc1];
end
