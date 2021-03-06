function emaxDAE()
% emaxDAE schedules PK-PD events, using a DAE solver. Drug effect should be
% observed within a given observable time.
% Infusion rate is a control variable.

% PK-PD parameters
p.Rate = 1; % infusion rate (mg/h)
p.V = 6; % volume of distribution (L)
p.k = 0.1; % elimination rate constant (1/h)
p.Emax = 20; % max effect (Hill equation)
p.E50 = 10; % drug concentration, linked to 50% mx effect (mg/l)

% call ODE

% speficy the Mass matrix
M = [1 0
     0 0];
    
options = odeset('Mass', M, 'RelTol',1e-6, 'AbsTol', 1e-6);
tspan = [0 25];

% initial value of the dependent variables
c0 = 0;
E0 = hillEffect(c0, 0, p.Emax, p.E50, 1);
[t,c] = ode15s(@derivalgeb, tspan, [c0 E0], options, p);

% generate graphs
set(axes, 'FontSize', 14)
plot(t, c(:,1), 'Color', 'k', 'LineWidth', 3)
hold on;
plot(t, c(:,2), 'Color', 'b', 'LineWidth', 3)
plot(t, rate(t,p),'Color', 'r', 'LineWidth',2)
grid minor;
grid on;
xlabel('Time [h]');
ylabel({'Drug Effect', 'Concentration [mg/L]',...
    'Infusion Rate [mg/h]'});
legend('conc', 'effect', 'rate')
% save graphs
print('-r600', '-dtiff', 'emaxDAE')

end

function massdydt = derivalgeb(t, y, p)
%massdyst = derivalgeb(t,y,p) calculates DAE model at points
%defined by the vector of dependent variables 'y', time 't', and
%with parameters 'p'. 

massdydt = [rate(t,p)/p.V - p.k*y(1);
    hillEffect(y(1), 0 , p.Emax, p.E50,1) - y(2)];
end

function r = rate(t, p)
% rate specify infusion rate within obervable time
%R = rate(t,p) calculates the infusion rate value at time 't'

tR = [0 1; 8 9; 15 16; 22 23]'; % rate times (begin end;...)
r = 0;
for ir = tR
    r = r + p.Rate * (heaviside(t-ir(1)) - heaviside(t-ir(2)));
end

end

function E = hillEffect(c, E0, Emax, EC50, n)
    E = E0 + Emax.*c.^n./(EC50.^n+c.^n);
end

