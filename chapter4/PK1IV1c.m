function PK1IV1c

% explicit solution for drug amounts and concentrations using
% one-compartment PK model.

clear all;
close all;
clc

p.V = 5.0;
p.k10 = 5.0;
dose = 100;
timeEnd = 10;

t = 0:0.1:timeEnd;
c = dose/p.V*exp(-p.k10*t);
a = c*p.V;
plot(t, a, '-r', t, c, 'b-', 'lineWidth', 2)
legend('Drug amounts [mg]', 'Drug concentrations [mg/L]')
xlabel('Time (h)')
ylabel('Value')
end