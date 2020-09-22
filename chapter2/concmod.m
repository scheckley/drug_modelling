function concmod(regimenType, numSubjects)

%concmod Model of Tumor Growth
% concmod(regimentype, numsubjects) models tumour growth and creates graphs
% of the tumour growth and epidermis over time. Regimentype must be either
% 'intermitent' or 'continuous'. Numsubjects is the number of subjects.

% example: model tumour growth with an intermittent regimen and 10
% subjects: concmod('intermitent',10)

% R. Gieschke and D. Serafin, Development of Innovative Drugs via Modeling
% with MATLAB, DOI: 10.1007/978-3-642-39765-3_2,
% Copyright Springer-Verlag Berlin Heidelberg 2014

% check input arguments and provide default values if necesary.

    narginchk(0,2);
    
    if nargin < 2
        numSubjects = 10;
    end
    
    if nargin < 1
        regimenType = 'intermittent';
    end
    
    % reset the random number generator to default
    rng default;
    
    % calculate dosing times and amount based on regimen   
    [doseTimes, doseAmount] = doseSchedule(regimenType);
    
    % set up figures for plotting
    tumorFigure = figure;
    tumorAxes = axes;
    xlabel('time [h]')
    ylabel('Numer of tumor cells (relative to baseline)')
    title('Tumor growth')
    xlim([0, doseTimes(end)]); hold on; grid on;
    epidermisFigure = figure;
    epidermisAxes = axes; set(epidermisAxes, 'FontSize', 14)
    xlabel('time [h]')
    ylabel('% change (relative to baseline)')
    title('Epidermis')
    xlim([0, doseTimes(end)]); hold on; grid on; ylim([50, 100]);
    
    %set([tumorFigure; epidermisFigure],...
    %    'units', 'normalized', ...
    %    {'Position'}, {[0.1, 0.5, 0.3, 0.4]; [0.6, 0.5, 0.3, 0.4]});
    
    
    % Simualte the system for each subject
    for subjectID = 1:numSubjects
        
        % initialize parameters and set initial conditions
        p = initializeParams;
        y0 = [doseAmount, p.c10, p.c20, p.n0, p.pc0, p.dc0, p.sc0];
        
        % pre-allocate variables to store the outputs
        timePoints = [];
        tumorGrowth = [];
        epidermis = [];
        
        % simulate the system for each treatment period
        for dose = 1:(length(doseTimes)-1)
            
            % set time interval for this treatment period
            tspan = [doseTimes(dose), doseTimes(dose+1)];
            
            % run the integrator
            [t,y] = ode45(@derivatives, tspan, y0, [], p);
            
            % record the outputs for plotting
            timePoints = [timePoints; t];
            tumorGrowth = [tumorGrowth; y(:,4)/p.n0];
            epidermis = [epidermis;100*y(:,7)/p.sc0];
            
            % reset initial conditions for the next treatment period
            % and add the next dose
            y0 = y(end,:);
            y0(1) = y0(1) + doseAmount;
        end
        
        % plot the results for this subject
        plot(tumorAxes, timePoints, tumorGrowth, 'Color', 'black')
        plot(epidermisAxes, timePoints, epidermis, 'Color', 'black')
        drawnow
    end
end