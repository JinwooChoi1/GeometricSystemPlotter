function pdot=odestepconstrainedoptimal(p,s,n,dimension,direction,nfparam)
    %%%%%%%%%%%%%
    % 
    % 
    %
    % Inputs:
    %
    % p: Matrix containing the Fourier series coefficients
    % s: System file which contains the connection vector field, CCF's and
    %   metric data
    % dimension: Indicates the number of shape variables of the system.
    % n: The number of points desired in a direct transcription parametrization
    %   of the gaits
    % costfunction: costfunction type to optimize for
    % lb: Lower bound of shape variables for each point which is obtained from the grid inside
    %   which an optimal gait is desired
    % ub: Upper bound of shape variables for each point which ll'l;is obtained from the grid inside
    %   which an optimal gait is desired
    % method: 'Hessian' - Based on the null vector of Hessian, the result
    %   vector is evaluated. Since the vector always point to the next step
    %   gait(stable), ode45 tend to skip some steps.
    %   'Perturb' - Based on the stable point of the gradient of
    %   Lagrangian function, the result vector is evaluated. At first, the
    %   vector does not point to the next step gait but it will fall into
    %   the stable point(next step) eventually. So, ode45 will try to
    %   calculate every step's solution.
    %
    % Outputs:
    %
    % pdot : the resulting vector is the negative gradient of displacement
    %   projected onto the null space of Hessian.
    %%%%%%%%%%%%%
    
    direction = [3 1];
    
    % Because of ODE solver, size of p is the number of fourier coefficient*dimension x 1
    p=reshape(p,[nfparam dimension]);

    disp1 = struct();
    disp2 = struct();
    stroke = struct();

    % when optimizing the steering gait, consider x and theta direction's
    % Hessian and gradient
    [jacobfourier1,temp_disp,temp_stroke] = evaluate_jacobian_fourier(p,s,n,dimension,direction(1));
    jacobfourier2 = evaluate_jacobian_fourier(p,s,n,dimension,direction(2));

    disp1.f = temp_disp(direction(1));
    disp2.f = temp_disp(direction(2));
    stroke.f = temp_stroke;

    % Calculate the Hessian of efficiency in the x direction.
    hessfourier1 = numerhessfourier(p,s,n,dimension,direction(1));
    hessfourier2 = numerhessfourier(p,s,n,dimension,direction(2));
    disp1.hf = cell2mat(hessfourier1.disp);
    disp2.hf = cell2mat(hessfourier2.disp);
    stroke.hf = cell2mat(hessfourier2.stroke);

    %% Record the fourier coefficient at step-optimal gaits.
    global currentDisp currentCost;

    currentCost = stroke.f;
    currentDisp = temp_disp;

    if strcmpi(s.costfunction,'torque') ||...
            strcmpi(s.costfunction,'covariant acceleration') ||...
            strcmpi(s.costfunction,'acceleration coord') ||...
            strcmpi(s.costfunction,'power quality')
        jacobfourier1.disp(end,:) = [];
        jacobfourier2.disp(end,:) = [];
        jacobfourier1.stroke(end,:) = [];
    end
    
    % Reshape jacobdisp and jacobstorke for ODE Solver.
    disp1.gf  = reshape(jacobfourier1.disp,[(nfparam-1)*dimension 1]);
    disp2.gf = reshape(jacobfourier2.disp,[(nfparam-1)*dimension 1]);
    stroke.gf = reshape(jacobfourier1.stroke,[(nfparam-1)*dimension 1]);

    lagr = evaluate_lagrangian(stroke,disp1,disp2);

    projgf = disp1.gf;
    optdir = -1;

    % Find the null space of Hessian and the null space of other constraint.
    U = null(lagr.hf,sqrt(0.1));
    V = null(disp2.gf.',sqrt(0.1));

    UV = [U -V];
    w = null(UV);

    nullvec = U*w(1:size(U,2),:);

    pdot = zeros(size(lagr.hf,1),1);

    % Project the gradient vector onto the null space
    for i = 1:size(nullvec,2)
        pdot=pdot+optdir*(nullvec(:,i).'*projgf)/norm(nullvec(:,i))^2*nullvec(:,i);
    end

    if isnan(pdot)
        pdot = zeros(size(pdot));
    end

    % The frequency should not change.
    pdot = reshape(pdot,[nfparam - 1 dimension]);
    pdot = [pdot; zeros(1,dimension)];
    pdot = reshape(pdot,[nfparam*dimension 1]);

    if(~isreal(pdot))
        pdot = real(pdot);
    end

end