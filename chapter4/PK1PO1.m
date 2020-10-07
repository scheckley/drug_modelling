function PK1PO1()

% solve one compartment model with first-order absorption and elimination.

p.k01 = 0.5; %1/h
p.V = 1.0; % L
p.CL = 2.0; % L/h
p.F = 1; % L/h
dose = 100 * p.F; % mg

timeEnd = 12;
[t,y] = ode45(@derivatives, 0:0.01:timeEnd, [dose; 0], [], p);

[AX, H1, H2] = plotyy(t,y(:,1), t,y(:,2), 'plot');
% set colours for plots
set(H1, 'LineWidth', 2, 'Color', 'blue');
set(H2, 'LineWidth', 2, 'Color', 'red');

% set colours for amount(t)
set(get(AX(1), 'Ylabel'),'String',...
    'Drug Amount at Site of Absorption', 'Color', 'blue',...
    'FontSize', 15)
set(AX(1), 'YColor', 'blue', 'FontSize', 15)
% set colours for concentration(t)
set(get(AX(2), 'Ylabel'), 'String', 'Drug Concentrations', 'Color',...
    'red', 'FontSize', 15)
set(AX(2), 'Ycolor', 'red', 'FontSize', 15)
title('One-compartment - Oral input')
xlabel('Time [h]')
legend('drug amount at size of absorption [mg]',...
    'drug concentrations [mg/L]')

end

function dydt = derivatives(~,y,p)
dy1 = -p.k01*y(1);
dy2 = p.k01*y(1) - p.CL/p.V*y(2);
dydt = [dy1; dy2];
end