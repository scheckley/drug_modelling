function pbpk1
% basic pbpk model representing blood and 1 organ.

p.Q = 2;
p.V = 5;
p.R = 1;
p.Vb = 10;
p.C = 0;
p.Cb = 10;

timeEnd = 12;
[t,y] = ode45(@derivatives, [0 timeEnd], [p.Cb; p.C], [], p);
plot(t,y)
hold on
p.R = 2;
[t,y] = ode45(@derivatives, [0 timeEnd], [p.Cb; p.C], [], p);
plot(t,y)

xlabel('Time [h]')
ylabel('Drug concentration [mg/L]');
legend('Blood', 'Organ', 'Location', 'South')
text(10, 7.2, 'R_{Organ} = 1')
text(10, 9.2, 'R_{Organ} = 2')
text(10, 4.2, 'R_{Organ} = 2')

end


function dydt = derivatives(~,y,p)

Cb = y(1);
C = y(2);
Ca = Cb;
Cv = C/p.R;
dydt(1) = p.Q * (Cv-Ca) / p.Vb;
dydt(2) = p.Q * (Ca-Cv) / p.V;
dydt = [dydt(1); dydt(2)];
end
