function E = hillEffect(c, E0, Emax, EC50, n)

    % hillEffect computes drug effect based on the Hill equation.
    % E = hillEffect(c, E0, Emax, EC50, n) computes drug effect E, based on
    % the Hill equation as a function of concentration c, and with with the
    % following concentration-response parameters:
    % E0 = baseline response, Emax = maximum effect, EC50 = concentration
    % related to 50% of max. effect, and n = Hill coefficient of
    % sigmoidicity.

    E = E0 + Emax.*c.^n./(EC50.^n+c.^n);
end