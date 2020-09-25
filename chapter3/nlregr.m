function nlregr()

% estimate PK parameters
% uses plot_ci from http://www.mathworks.com/matlabcentral/fileexchange/31752-plotci
% under BSD license. Copyright Zbigniew (2011).

observations = [0 1 2 3 4 5 6 12 %time points
                0 1 1.5 1.6 1.5 1.4 1.3 1.4]'; %drug amount

tObs = observations(:,1);
aObs = observations(:,2);

% beta estimation
options = statset('Robust', 'on', 'WgtFun', 'bisquare');
[beta, resid, J] = nlinfit(tObs, aObs, @model, [1 2], options);

%objective function
w = (abs(resid)<1).*(1-resid.^2).^2; %weight function 'bisquare'
ofv = resid' * diag(w) * resid;

% confidence interval for model parameters
paramCI = nlparci(beta, resid, 'jacobian', J);

disp([beta,'paramCI']);

% confidence interval for predictions
tPred = 0:0.02:12;

[aPred, delta] = nlpredci(@model, tPred, beta, resid, 'jacobian',...
    J, 'simopt', 'on');

[~, iInd] = ismember(tObs, tPred);

% show predictions
disp(['time', 'obs','Pred']);
disp(num2str([tObs aObs, aPred(iInd)]))
predCI = [aPred aPred+delta aPred-delta];
set(axes, 'FontSize', 15)
plot(tObs, aObs, 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on
hold on
rgbfill = [0.3047 0.6992 0.9102];
rgbline = [0.9414 0.9414 0.3359];

% call plot_ci
plot_ci(tPred, predCI, 'PatchColor', rgbfill, 'PatchAlpha', 0.5,...
    'MainLineWidth', 2, 'MainLineColor', rgbline, 'LineWidth', 1,...
    'LineStyle','-', 'LineColor', 'k');

legend('Observation', 'Prediction CI', 'Regression line',...
    'Location', 'Best');

xlabel('Time')
ylabel('amount')
title(['Fitting:','ka=', num2str(beta(1)), ' k=',...
    num2str(beta(2)),' OFV(least squares)=', num2str(ofv)])
print('-r1000', '-dtiff', 'nlregression');
end

function dadt = derivatives(~, a, p)
% compute ode's

% PK one-compartment model (SC administration)
dadt = [-p.ka * a(1);               %drug amount - depot
         p.ka * a(1) - p.k * a(2)]; % drug amount - compartment
end

function a = model(beta, t)
%model prediction of the drug amount
% A = model(beta, t) computes the prediction for the drug amount 'A'
% in one-compartmetn model at time points given as vector 'T',
% and with model parameters 'beta'

dose = 100;
amount0 = 0;
tspan = [0 12];
a0 = [dose amount0];
p.ka = beta(1); % absorption rate
p.k = beta(2); % elimination rate
sol = ode45(@derivatives, tspan, a0, [], p);
a12 = deval(sol, t)';
a = a12(:,2);
end




