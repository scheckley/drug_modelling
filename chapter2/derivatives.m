function dydt = derivatives( ~, y, p)

% Compute the right hand side of the ODE
% dydt = derivatives(t,y,p)

    amount = y(1); %drug amount
    c1 = y(2); % parent drug concentration
    c2 = y(3); % active metabolite concentration
    n = y(4); % tumour growth
    pc = y(5); % cells in the proliferative compartment
    dc = y(6); % cells in the differentiated compartment
    sc = y(7); % cells in the stratum corneum compartment
    
    % pk model
    dAmountdt = -p.k01*amount;
    dc1dt = p.k01/p.V1*amount - (p.CL12/p.V1 + p.CL10/p.V1)*c1;
    dc2dt = p.CL12/p.V2*c1 - p.CL20/p.V2*c2;
    
    % efficacy model
    kkill = hillEffect(p.tumFactor*c2, 0, p.kkillMax, p.kkill50, 1);
    dndt = (p.k0*log(p.nn00/n) - (p.k + kkill))*n;
    
    % toxicity model
    DS = hillEffect(p.skFactor*c2, 0, p.DSMax, p.DS50,1);
    br1 = p.br0 * (1-DS); % birth rate
    dpcdt = br1 - p.kpc*pc;
    ddcdt = p.kpc*pc - p.kdc*dc;
    dscdt = p.kdc*dc - p.ksc*sc;
    
    % Derivatives vector of ODE system
    
    dydt = [dAmountdt; %drug amount
            dc1dt; % parent drug concentration
            dc2dt; % active metabolite concentration
            dndt; % tumour growth
            dpcdt; % changes in the proliferative compartment
            ddcdt; % changes in the differentiated compartment
            dscdt;  % changes in stratum corneum compartment
            ];
end