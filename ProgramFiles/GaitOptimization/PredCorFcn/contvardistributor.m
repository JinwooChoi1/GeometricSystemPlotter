function [p,sigma,lambda,mu,c] = contvardistributor(y,gopt)
    %%%%%%%%%%%%%
    % This function distribute continuation parameters to gait parameters,
    % Lagrange multipliers, and continuation variables.
    %
    % Inputs:
    %
    % y: The continuation variables. It contains [p; sigma; lambda; mu; c].
    %
    % Outputs:
    % p: Matrix containing the Fourier series coefficients(Gait parameter).
    % sigma, lambda, mu: The multipliers for the cost, equality constraint,
    %   and the inequality constraint, respectively.
    % c: The scalar parameter for continuation. In this problem, it is the 
    %   step size of the gait.
    %%%%%%%%%%%%%
    nptotal = gopt.nftotal - gopt.dimension;
    
    if isfield(gopt,'neq')
        neq = gopt.neq;
    else
        neq = 1;
    end

    p = y(1:nptotal);
    sigma = y(nptotal+1);
    lambda = y(nptotal+2:nptotal+neq+1);
    mu = y(nptotal+2+neq:end-1);
    c = y(end);
end