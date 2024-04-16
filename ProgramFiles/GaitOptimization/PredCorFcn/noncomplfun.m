function q = noncomplfun(q,gopt)

    dimension = gopt.dimension;
    nf = gopt.nfparam;
    A = gopt.A;
    b = gopt.b;

    [p,sigma,lambda,mu,c] = contvardistributor(q,gopt);

    if ~isempty(mu)
        A(:,nf*(1:dimension)) = [];
    
        mu(A*p-b < -1e-3) = 0;
        mu(abs(A*p-b) < 1e-3) = 1e-3;
        mu(mu < 0) = 1e-3;
    end
    q = [p; sigma; lambda; mu; c];
end