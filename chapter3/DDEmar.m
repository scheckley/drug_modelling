function DDEmar
% DDEmar solves a PK model with time delay
% DDEmar is a simple example that solves a PK model
% described as a DDE. The simulation is repeated
% for several delay values with comparison plot.

% PK parameters
p.F = 1;
p.CL = 2.16;
p.V = 12.6;
p.ka = 0.1;

% parameters for calling dde23
p.c0 = 0;
p.A0 = 100;
tspan = [0 60];
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-6);
alag = {3 2 1}; % possible delays
a = axes;
set(a, 'FontSize', 15)
for si = struct('alag', alag, 'col', {'b', 'g', 'r'})
    sol = dde23(@derivatives, si.alag, @history, tspan, options, p);
    t = sol.x;
    y = sol.y;
    plot(t, (y(2,:)), si.col, 'LineWidth', 2);
    hold on
end

% compare with analytical solution without delay (alag=0)
Dose = p.A0;
F = p.F;
ke = p.CL/p.V;
ka = p.ka;
V = p.V;
y_ana = (F*Dose/V)*(ka/(ka-ke))*(exp(-ke*t) - exp(-ka*t));
plot(t, y_ana, 'k', 'LineWidth', 2);

legend(['t_{lag}=' num2str(alag{1})], ['t_{lag}='...
    num2str(alag{2})], ['t_{lag}=' num2str(alag{3})], ...
    't_{lag}=0', 'Location', 'Best')

title('Concentration-Time Course');
xlabel('Time [h]')
ylabel('Concentration [\mug/L]')
print('-r900', '-dtiff', 'ddemar')

end

function dydt = derivatives(~, y, Z, p)

% compute dydt. Matrix Z corresponds to all possible delays
% and contains solutions 'y' at the delayed argument. 
% Thus, 'Z' has a number of rows equal to the length of 'y',
% and the number of columns equal to the number of delays.
% for example, Z(2,1) means y(2) at time t-clag where clag is a delay.

% PK model
a = y(1); % amount
%c = y(2) %y(2) is not used directly
clag = Z(2,1);
dadt = -p.ka*a;
dcdt = (p.ka/p.V*p.F*a - p.CL/p.V*clag);
dydt = [dadt
    dcdt];

end

function yhist = history(~,p)
% history provides the history of the solution
% yhist = history(t, p) specifies the time profiles of dependent
% variables at times 't' prior to the initial point. 
% This function is a form of 'initial condition' for the DDE.
% 'p' stands for additional parameters which can be passed to the function.

yhist = [p.A0, p.c0];
end




