function params = initializeParams

% initalizeParams creates initial values for the model parameters
% params = initializeParams returns a structure 'params' containing initial
% values for model parameters, including variability when necessary.

    % PK parameters
    pkCV = 0.3;
    lim1 = 0.3;
    lim2 = 5.0;
    
    k01 = 0.7; % 1st order elimination
    params.k01 = variability('varnorm', k01, k01 * pkCV, 1, lim1);
    CL12 = 10; % metabolite clearance
    params.CL12 = variability('varnorm', CL12, CL12 * pkCV, 1, lim2);
    CL10 = 80; % parent drug clearances
    params.CL10 = variability('varnorm', CL10, CL10 * pkCV, 1, lim2);
    CL20 = 60; % metabolite elimination clearances
    params.CL20 = variability('varnorm', CL20, CL20 * pkCV, 1, lim2);
    V1 = 30; % metabolite volumes of distribution
    params.V1 = variability('varnorm', V1, V1 * V1, 1, lim2);
    V2 = 150; % metabolite elimination volumes of distribution
    params.V2 = variability('varnorm', V2, V2 * pkCV, 1, lim2);
    params.c10 = 0; % initial value: parent drug concentration
    params.c20 = 0; % initial value: active metabolite concentration
    
    % efficacy parameters (tumor growth parametes)
    params.k0 = 4.2e-5; % doubling time 105 days
    params.k = 0; % natural elimination rate
    params.tumFactor = 12; % tumour factor
    params.nn00 = 1e12; % tumour growth limiting value
    params.n0 = 1e9; % tumour growth initial value
    params.kkillMax = 0.05; % max effect (hill equation)
    kkill50 = 100; % concentration linked to 50% of maximum effect
    kkill50CV = 1.5;
    params.kkill50 = variability('varlog',...
        kkill50, kkill50CV*kkill50, 1);
        
    % Toxicity parameters (epidermis)
    params.skFactor = 8.0; % toxicity factor
    params.DSMax = 0.8; % maximum effect (hill equation)
    DS50 = 10.0; % concentration linked to 50% of maximum effect
    DS50CV = 0.8; 
    params.DS50 = variability('varlog',...
        DS50, DS50CV * DS50, 1);
        
    % proliferative compartment (pc)
    ttpcd = 6; % cell cycle time in pc (days)
    ttpc = ttpcd * 24; % cell cycle time in pc (hours)
    params.kpc = 1/ttpc; % app. elimination rate at steady state
    growthFactor = 1.0; % growth factor in pc
    pcTot = 27000; % number of cells in pc
    params.pc0 = growthFactor*pcTot; % cells in pc at time 0
    params.br0 = params.pc0/ttpc; % initial birth rate
    
    % differentiated compartment (dc)
    ttdcd = 9; % transit time in dc (days)
    ttdc = ttdcd * 24; % transit time in dc (hours)
    params.kdc = 1/ttdc; % app. elimination rate at steady state
    params.dc0 = params.br0*ttdc; % cells in dc at time 0
    
    % stratum corneum compartment (sc)
    ttscd = 7; % transit time in sc (days)
    ttsc = ttscd*24; % transit time in sc (hours)
    params.ksc = 1/ttsc; % app. elimination rate at steady state
    params.sc0 = params.br0*ttsc; % cells in sc at time 0
    
end