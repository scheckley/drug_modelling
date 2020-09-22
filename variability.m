function x = variability(distribution, m, s, num, lim)
    % variability: generate variability
    % returns a vector of random nubmers for varaibility, specified by
    % 'distribution' with the mean value 'm' and standard deviation 's'.
    % 'num' is the length of the vector 'x'. 'distribution' can be 'varlog'
    % or 'varnorm': 'varlog' generates random numbers x-logN(n,s). m
    % stands for mean value of the log normal distribution. 'varnorm'
    % generates normally distributed random numbers 'x', where x ~ n(m,s),
    % right-censored by 'lim'
    
    switch distribution
        case 'varlog'
            %lognormal distribution
            mu = log(m^2/sqrt(m^2 + s^2));
            sigma = sqrt(log(1+(s/m)^2));
            x = lognrnd(mu, sigma, 1, num);
        case 'varnorm'
            %lognormal distribution
            mu = m;
            sigma = s;
            x = max(normrnd(mu, sigma, 1, num), lim);
        otherwise
            % placeholder for any other distribution
    end
    
end    