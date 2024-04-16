function pdot=odestepoptimal(p,s,n,dimension,direction,nf)%,lb,ub,writerObj)
    %%%%%%%%%%%%%
    % This function calculates the relationship vector between one optimal
    % gait and next optimal gait with respect to the coefficients obtained
    % by the fourier series parametrization.
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
    %
    % Outputs:
    %
    % pdot : the relationship vector between one optimal gait and next
    % optimal gait.
    %%%%%%%%%%%%%

    % Because of ODE solver, size of p is the number of fourier coefficient*dimension x 1
    p=reshape(p,[nf-1 dimension]);
    p = [p; 2*pi*ones(1,dimension)];

    nftotal = (nf-1)*dimension;

    disp = struct();
    stroke = struct();

    % Derive Hessian and gradient
    [jacobfourier,temp_disp,temp_stroke] = evaluate_jacobian_fourier(p,s,n,dimension,direction);
    disp.f = temp_disp(direction);
    stroke.f = temp_stroke;

    hessfourier = numerhessfourier(p,s,n,dimension,direction);
    disp.hf = cell2mat(hessfourier.disp);
    stroke.hf = cell2mat(hessfourier.stroke);

    %% Record the fourier coefficient at step-optimal gaits.
    global currentDisp currentCost;

    currentCost = stroke.f;
    currentDisp = disp.f;

    %% calcuating pdot so that lagrange equation is zero.
    if strcmpi(s.costfunction,'torque') ||...
            strcmpi(s.costfunction,'covariant acceleration') ||...
            strcmpi(s.costfunction,'acceleration coord') ||...
            strcmpi(s.costfunction,'power quality')
        jacobfourier.disp(end,:) = [];
        jacobfourier.stroke(end,:) = [];
        jacobfourier.eff(end,:) = [];
    end

    % Reshape jacobdisp and jacobstorke for ODE Solver.
    disp.gf = reshape(jacobfourier.disp,[nftotal 1]);
    stroke.gf = reshape(jacobfourier.stroke,[nftotal 1]);

    lagr = evaluate_lagrangian(stroke,disp);
    
    projgf = disp.gf;
    optdir = 1;

    % Find the null space of Hessian. Its tolerance is 0.1.
    nullvec = null(lagr.hf,0.2);
    pdot = zeros(nftotal,1);

    % Project the gradient vector onto the null space
    for i = 1:size(nullvec,2)
        pdot=pdot+optdir*(nullvec(:,i).'*projgf)/norm(nullvec(:,i))^2*nullvec(:,i);
    end

    if isnan(pdot)
        pdot = zeros(size(pdot));
        warning("The solution contains NaN.");
    end
end