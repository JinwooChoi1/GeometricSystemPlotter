function [H, JH, projvec, stop] = steeringcontfunc(y,s,gopt)%,lb,ub,writerObj)
    %%%%%%%%%%%%%
    % This function calculates the relationship vector between one optimal
    % gait and next optimal gait with respect to the coefficients obtained
    % by the fourier series parametrization.
    %
    %
    % Inputs:
    %
    % y: The continuation variables. It contains [p; sigma; lambda; mu; c].
    %   p : the fourier descriptor of the gait.
    %   sigma,lambda,mu : the lagrange multiplier
    %   c : the continuation variables.
    % s: System file which contains the connection vector field, CCF's and
    %   metric data
    % gopt: struct generated by optimalgaitoptions().
    %
    % Outputs:
    %
    % H : The continuation function.
    % JH : Jacobian of the continuation function.
    % projvec : The predictor guidance vector.
    % stop : stop condtion. If c is small enough, stop the pc iteration.
    %
    %%%%%%%%%%%%%

    dimension = gopt.dimension;
    np = gopt.nfparam;
    anchor = gopt.anchor;

    nptotal = (np-1)*dimension;    

    if isfield(gopt,"neq")
        neq = gopt.neq;
        c02 = gopt.c02;
    else
        neq = 1;
    end

    p = contvardistributor(y,gopt);

    p = reshape(p,[np-1 dimension]);
    p = [p; 2*pi*ones(1,dimension)];
    disp1 = struct();
    disp2 = struct();
    stroke = struct();

    % when optimizing the steering gait, consider x and theta direction's
    % Hessian and gradient
    gopt.issubopt = false;
    [jacobfourier1,~,temp_stroke,disp1.f] = evaluate_jacobian_fourier(p,s,gopt);
    hessfourier1 = numerhessfourier(p,s,gopt);
    gopt.issubopt = true;
    [jacobfourier2,~,~,disp2.f] = evaluate_jacobian_fourier(p,s,gopt);
    hessfourier2 = numerhessfourier(p,s,gopt);
    gopt.issubopt = false;
    stroke.f = temp_stroke;

    % Reshape jacobdisp and jacobstorke for ODE Solver.
    disp1.gf = reshape(jacobfourier1.disp,[], 1);
    disp2.gf  = reshape(jacobfourier2.disp,[], 1);
    stroke.gf = reshape(jacobfourier1.stroke,[], 1);

    % Calculate the Hessian of efficiency in the x direction.
    disp1.hf = cell2mat(hessfourier1.disp);
    disp2.hf = cell2mat(hessfourier2.disp);
    stroke.hf = cell2mat(hessfourier1.stroke);

    eff1 = fraction_derivative(disp1,stroke);
    eff2 = fraction_derivative(disp2,stroke);

    [eff1, eff2] = normalize_eff(eff1, eff2, gopt.anchor);
    
    if neq == 2
        obj = deal(stroke);
        [const{2},const{1}] = atan2_derivative(eff2,eff1);
        const{1}.f = const{1}.f-c02;
%         projvec = [zeros(length(y)-1,1);-1];
        projvec = [-stroke.gf; zeros(length(y)-nptotal-1,1); -1];

        if const{2}.f < 0.1
            stop = 1;
        else
            stop = 0;
        end

    else
        [obj,const] = atan2_derivative(eff2,eff1);
        
        projvec = [const.gf; zeros(length(y)-nptotal-1,1);1];

        if (gopt.contdir == -1)
            projvec = -projvec;
            stop = (const.f < 0.1);
        else
            stop = (const.f > pi/2);
        end

    end
    
    [H,JH] = constructHfunc(y,obj,const,gopt);
end