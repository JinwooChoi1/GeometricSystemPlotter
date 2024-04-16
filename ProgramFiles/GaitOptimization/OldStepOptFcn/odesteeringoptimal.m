function pdot=odesteeringoptimal(p,s,n,dimension,direction,nf)
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
    p = reshape(p,[nf-1 dimension]);
    p = [p; 2*pi*ones(1,dimension)];

    nftotal = (nf-1)*dimension;

    disp = struct();
    stroke = struct();    
    
    % when optimizing the steering gait, consider x and theta direction's
    % Hessian and gradient
    [jacobfourierx,temp_disp,temp_stroke,dispx.f] = evaluate_jacobian_fourier(p,s,n,dimension,direction(1));
    stroke.f = temp_stroke;
    [jacobfourier,~,~,disp.f] = evaluate_jacobian_fourier(p,s,n,dimension,direction(2));
    % Calculate the Hessian of efficiency in the x direction.
    hessfourierx = numelhessfourier(p,s,n,dimension,direction(1));
    hessfourier = numelhessfourier(p,s,n,dimension,direction(2));
    dispx.hf = cell2mat(hessfourierx.disp);
    disp.hf = cell2mat(hessfourier.disp);
    stroke.hf = cell2mat(hessfourier.stroke);
    
    % Reshape jacobdisp and jacobstorke for ODE Solver.
    disp.gf = reshape(jacobfourier.disp,[nftotal 1]);
    stroke.gf = reshape(jacobfourier.stroke,[nftotal 1]);
    dispx.gf  = reshape(jacobfourierx.disp,[nftotal 1]);

    if (rotOptMod == 0)
        obj = fraction_derivative(dispx,stroke);
        const = fraction_derivative(disp,stroke);
        optdir = 1;
    else
        obj = fraction_derivative(disp,stroke);
        const = fraction_derivative(dispx,stroke);
        optdir = -1;
    end

    lagr = evaluate_lagrangian(obj,const);

    % Find the null space of Hessian. Its tolerance is 0.2

    nullvec = null(lagr.hf,0.2);

    pdot = zeros(size(lagr.hf,1),1);

    % Project the gradient vector onto the null space
    for i = 1:size(nullvec,2)
        pdot=pdot+optdir*(nullvec(:,i).'*const.gf)/norm(nullvec(:,i))^2*nullvec(:,i);
    end

    if isnan(pdot)
        pdot = zeros(size(pdot));
        warning("The solution contains NaN.");
    end
end