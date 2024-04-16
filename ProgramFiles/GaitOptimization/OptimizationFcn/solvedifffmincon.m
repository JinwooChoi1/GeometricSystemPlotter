function [f,g]=solvedifffmincon(y,s,gopt)%,lb,ub,writerObj)
    %%%%%%%%%%%%%
    % This function calculates efficiency (or displacement, if
    % that is the objective function) and its gradient with respect to the coefficients obtained
    % by the fourier series parametrization
    %
    % Inputs:
    %
    % y: Matrix containing the Fourier series coefficients
    % s: System file which contains the connection vector field, CCF's and
    %   metric data
    % dimension: Indicates the number of shape variables of the system
    % n: The number of points desired in a direct transcription parametrization
    %   of the gaits
    % direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
    % costfunction: costfunction type to optimize for
    %
    % Outputs:
    %
    % f: Objective function value (This is negative of efficiency by default, can be
    %   changed to displacement)
    % g: Gradient of the objective function
    %%%%%%%%%%%%%
    global bestCost bestDisp bestEff;

    dimension = gopt.dimension;
    if (gopt.issubopt)
        direction = gopt.subdirection;
    else
        direction = gopt.direction;
    end

    if nargout > 1
        [jacobfourier,avgd_vel,totalstroke] = evaluate_jacobian_fourier(y,s,gopt);
        lineint = avgd_vel(direction);

        if lineint < 0
            lineint = -lineint;
            jacobfourier.disp = -jacobfourier.disp;
            jacobfourier.eff = jacobfourier.disp/totalstroke-...
                lineint*jacobfourier.stroke/totalstroke^2;
        end
    else 
        s.costfunction = gopt.costfunction;
        [~, net_disp_opt, totalstroke] = evaluate_displacement_and_cost1(s,y);
        net_disp_opt = se2log(net_disp_opt);
        lineint = net_disp_opt(direction);

        if lineint < 0
            lineint = -lineint;
        end
    end


%% minimizing negative of efficiency(or displacement)
f=-lineint/(totalstroke); % Optimizing for displacement over cost
%     f = -lineint; % Optimizing for displacement only

if abs(f) > bestEff
    bestEff = abs(f);
    bestDisp = abs(lineint);
    bestCost = abs(totalstroke);
end

if nargout > 1
    % Return the gradient of efficiency plus row of zeros for frequency
    % terms for drag systems
    g= [-jacobfourier.eff;zeros(1,dimension)]; % Optimizing for displacement over cost
    % g = [-jacobfourier.disp;zeros(1,dimension)]; % Optimizing for displacement only
end


end