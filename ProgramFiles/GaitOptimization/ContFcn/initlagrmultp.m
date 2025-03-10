function [multp,c,c2] = initlagrmultp(s,y,mu0,gopt)

    dimension = gopt.dimension;
    direction = gopt.direction;
    nf = gopt.nfparam;
    A = gopt.A;
    b = gopt.b;

    if ~isfield(gopt,"neq")
        neq = 1;
    else
        neq = gopt.neq;
    end

    % Calculate the gradient of constraint and objective function.
    % Baby step optimization
    if ~isfield(gopt,'subdirection')
        [jacobfourier,avgd_vel] = evaluate_jacobian_fourier(y,s,gopt);

        obj.gf = reshape([jacobfourier.stroke; zeros(1,dimension)],[nf*dimension 1]);
        const.gf = [reshape([jacobfourier.disp; zeros(1,dimension)],[nf*dimension 1])];
        c = avgd_vel(direction);

    % Steering optimization
    else
        anchor = gopt.anchor;

        gopt.issubopt = false;
        [jacobfourier1,~,stroke.f,disp1.f] = evaluate_jacobian_fourier(y,s,gopt);
        gopt.issubopt = true;
        [jacobfourier2,~,stroke.f,disp2.f] = evaluate_jacobian_fourier(y,s,gopt);

        stroke.gf = reshape([jacobfourier1.stroke; zeros(1,dimension)],[nf*dimension 1]);
        disp1.gf = reshape([jacobfourier1.disp; zeros(1,dimension)],[nf*dimension 1]);
        disp2.gf = reshape([jacobfourier2.disp; zeros(1,dimension)],[nf*dimension 1]);

        eff1 = fraction_derivative(disp1,stroke);
        eff2 = fraction_derivative(disp2,stroke);  

        eff1.f = (eff1.f-anchor(2,1))/(anchor(1,1)-anchor(2,1));
        eff1.gf = eff1.gf/(anchor(1,1)-anchor(2,1));
        eff2.f = (eff2.f-anchor(1,2))/(anchor(2,2)-anchor(1,2));
        eff2.gf = eff2.gf/(anchor(2,2)-anchor(1,2));

        [obj,const] = atan2_derivative(eff2,eff1);

        if neq == 2
            c = obj.f;
            c2 = const.f;
        elseif neq == 1
            c = const.f;
            c2 = [];
        else
            error("The number of equality constraint is wrong.")
        end
    end

    % If it is the steering problem, it is a maximization problem.
    % If it is the baby step problem, it is a minimization problem. (KKT conditions)
    if (isfield(gopt,'subdirection')) && (neq == 1)
        signobj = -1;
    else
        signobj = 1;
    end

    % include the inequality multiplier if the inequality constraint is
    % violated or the optimization problem is the steering problem.
    if any(A*reshape(y,[], 1) - b >= -1e-2) || isfield(gopt,'subdirection')
        ineq.f = A*reshape(y,[], 1)-b;
        ineq.gf = A.';        
        
        % Set the number of equality constraint and inequality constraint.
        m1 = size(const.gf,2);
        m2 = size(ineq.gf,2);

        % The inequality multiplier should be positive.
        lA = blkdiag(zeros(m1),-eye(m2));
        lb = zeros(m1+m2,1);

        % Complemantary slackness.
        lAeq = diag([zeros(m1,1); (ineq.f < -1e-6)]);
        lbeq = zeros(m1+m2,1); 
        gopt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','none');
        if ~exist('mu0','var') || isempty(mu0)
            mu0 = zero(length(b),1);
        end
        lambda0 = [1; mu0];

        if neq == 2
            C = [const.gf obj.gf ineq.gf];
            d = signobj*stroke.gf;
            lA = blkdiag(0,lA);
            lb = [0; lb];
            lAeq = blkdiag(0,lAeq);
            lbeq = [0; lbeq];
            lambda0 = [1; lambda0];
        elseif neq == 1
            C = [const.gf ineq.gf];
            d = signobj*obj.gf;
        else
            error("The number of equality constraint is wrong.")
        end

        % Solve the duality problem. It can be expressed as a nonlinear 
        % least sqaure problem. However, use the quadratic programming for
        % computation speed (||Cx+d||^2).
        lambda = quadprog(C.'*C,C.'*d,lA,lb,lAeq,lbeq,[],[],lambda0,gopt);

        % Add the objective multiplier.
        multp = [signobj; lambda];
        multp = multp/norm(multp);
    else
        % Solve the duality problem. It can be expressed as a nonlinear 
        % least sqaure problem. The solution is simlply the moore-penrose
        % inverse.
        if neq == 1
            lambda = -pinv(const.gf)*obj.gf;
        else
            lambda = -pinv([const.gf obj.gf])*stroke.gf;
        end
        % Add the objective multiplier.
        multp = [1; lambda];
        multp = multp/norm(multp);
    end

end