function infus1inter()

    %infus1inter() one-compartment model with infuction. 
    % computes plasma concentration time course and
    % plots output. function has 2 parameters, volume (V) and clearance (CL)

    % PK parameter(s) specification
    p.V = 12.6; % volume of distribution
    p.CL = 2.16; % clearance

    
    % set solver options
    options = odeset;
    options.RelTol = 1e-6;
    options.AbsTol = 1e-6;
    
    % or just use ode23!
    
    % run the simulator
    
    
    [t,c] = ode45(@derivatives, 0:0.01:15, 0, options, p);

    % plot output

    subplot(3,1,1);
    plot(t, r(t), 'Color', 'k', 'LineWidth', 2)
    title(['One compartment model with infusion', 10], 'FontSize', 15)
    ylim([0 110])
    ylabel({'Infusion' 'Rate [mg/h]'}, 'FontSize', 12)
    subplot(3,1, 2:3);
    plot(t,c, 'Color', 'k', 'LineWidth',2)
    xlabel('Time [h]', 'FontSize',12)
    ylabel('Concentration [mg/L]', 'FontSize', 12)
    %save the plot
    print('-r900','-dtiff','infus1inter')
end

function dcdt = derivatives(t, c, p)
    dcdt = r(t)/p.V - p.CL/p.V*c;
end


function rt = r(t)
    % R specifies infusion rate
    % RT = R(t) returns infusion rate RT and time T

    rt = 200/2 * (t < 2);
end

