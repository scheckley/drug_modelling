function pbpk3
% PBPK model with 2 organs and drug elimination

% parameters
% flows
p.Q1 = 2;
p.Q2 = 1;
% volumes
p.Vb = 5;
p.V1 = 4;
p.V2 = 4;
% partition coefficients
p.R1 = 0.5;
p.R2 = 2;
% elimination parameters
p.CL = 1;
% initial values
p.cbi = 10;
p.c1i = 0;
p.c2i = 0;

timeEnd = 24;
[t, y] = ode45(@derivatives, [0 timeEnd], [p.cbi; p.c1i; p.c2i], [], p);
plot(t,y,'LineWidth', 2)
hold on
p.CL=0;
[t, y] = ode45(@derivatives, [0 timeEnd], [p.cbi, p.c1i, p.c2i], [], p);
plot(t,y,'-')
xlabel('Time [h]');
ylabel('Drug Concentrations [mg/L]');
legend('Blood','Organ 1 - non-eliminating',...
    'Organ 2 - eliminating')

end

function dydt = derivatives(~,y,p)

cb = y(1);
c1 = y(2);
c2 = y(3);
Q = p.Q1+p.Q2;
c1v = c1/p.R1;
c2v = c2/p.R2;
ca = cb;
cv = (c1v * p.Q1 + c2v * p.Q2) / Q;
rout = p.CL * c2;
dcbdt = Q * (cv - cb) / p.Vb;
dc1dt = p.Q1 * (ca - c1v) / p.V1;
dc2dt = (p.Q2 * (ca - c2v) - rout) / p.V2;
dydt = [dcbdt; dc1dt; dc2dt];

end