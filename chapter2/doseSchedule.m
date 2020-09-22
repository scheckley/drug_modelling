function [doseTimes, doseAmount] = doseSchedule(regimenType)

% doseSchedule creates a vector of dosing times and amounts of drug.
% [doseTimes, doseAmount] = doseSchedule(regimenType) returns a vector
% doseTimes, of dose times and a scalar doseAmount of the dosing amount.
% RegimenType must be either 'intermittent' or 'continuous'.

    treatmentWeeks = 12;
    cycleWeeks = 3;
    daysInWeek = 7;
    hoursInDay = 24;
    initialTime = 0;
    
    % treatment time [hours]
    endTime = treatmentWeeks * daysInWeek * hoursInDay;
    
    % cycle time [hours]
    cycleTime = cycleWeeks * daysInWeek * hoursInDay;
    
    switch regimenType
        case 'intermittent'
            % intermittent treatment
            dailyDoses = 2; %BID
            doseAmount = 1500; % 1 dose in [mg]
            timeOnDrug = cycleTime - daysInWeek * hoursInDay;
        case 'continuous'
            % contunuous (over 12 weeks) treatment
            dailyDoses = 2; % BID
            doseAmount = 1000; % 1 dose in [mg]
            timeOnDrug = cycleTime - hoursInDay/dailyDoses;
        otherwise
            % case 'other' placeholder for other regimen(s)
            messageID = 'concmod:UnknownRegimen';
            messageStr = ['Unknown regimen ''%s''. The regimen must ',...
                'be either ''intermittent'' or ''continuous''.'];
            error(messageID, messageStr, regimenType);
            
    end
    
    % create vector of treatment times
    initialCycleTimes = ...
        initialTime : cycleTime : (endTime - cycleTime);
   
    doseTimesWithinCycle = ...
        initialTime : hoursInDay/dailyDoses : timeOnDrug;
    
    doseTimes = reshape(...
        repmat(doseTimesWithinCycle', 1, length(initialCycleTimes)) + ...
        repmat(initialCycleTimes, length(doseTimesWithinCycle), 1), ...
        1, length(initialCycleTimes) * length(doseTimesWithinCycle));
    doseTimes = [doseTimes, endTime];
    
end