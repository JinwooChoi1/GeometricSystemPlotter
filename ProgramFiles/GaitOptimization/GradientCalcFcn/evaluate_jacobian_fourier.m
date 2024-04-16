function [jacobfourier,avgd_vel,totalstroke,avgd_vel_dir] = evaluate_jacobian_fourier(f,s,gopt)
    %%%%%%%%%%%%%
    % This function calculates the gradient of displacement and cost function
    % with respect to the coefficients obtained by
    % the fourier series parametrization. It call evalute_jacobian.
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
    % jacobfourier: the struct which contains the gradient of functions.
    %   (jacobfourier.disp, jacobfourier.stroke,
    %   jacobfourier.eqi, jacobfourier.eff)
    % lineint : the displacement for the gait
    % totalstroke : the cost for the gait
    %%%%%%%%%%%%%

    n = gopt.npoints;
    dimension = gopt.dimension;
    
    %% Calculate the jacobian with respect to the discretized point on the gait cycle
    [jacobian,avgd_vel,totalstroke,avgd_vel_dir] = evaluate_jacobian(f,s,gopt);

    %% Check the system type.
    % for inertia-dominated systems, systemtype = 1.
    % for drag-dominated systems, systemtype = 0.
    if strcmpi(gopt.costfunction,'torque') ||...
            strcmpi(gopt.costfunction,'covariant acceleration') ||...
            strcmpi(gopt.costfunction,'acceleration coord') ||...
            strcmpi(gopt.costfunction,'power quality')
        systemtype = 1;
    else
        systemtype = 0;
    end

    %% chy tells us how much each point moves when a fourier series variable is changed
    chy=chy_generator(f,n,dimension);

    %% properly ordering gradients depending on wether lineint was negative or positive
    % if invert
    %     jacobian.disp=flip(jacobian.disp);
    %     jacobian.eqi=flip(jacobian.eqi);
    %     jacobian.stroke=flip(jacobian.stroke);
    %     jacobian.repuls=flip(jacobian.repuls);
    % end



    %% fourier series version of all gradients

    % We then obtain gradients in a fourier series parametrization by
    % projecting the gradients from the direct transcription space onto the
    % fourier coefficient space
    jacobfourier=struct();
    jacobfourier.disp = zeros(length(f)-1,dimension);
    jacobfourier.stroke = zeros(length(f)-1,dimension);
    jacobfourier.eqi = zeros(length(f)-1,dimension);

    for i=1:1:dimension
        for j=1:1:length(f)-1
            jacobfourier.disp(j,i)=chy{i}(j,:)*jacobian.disp(:,i);
            if systemtype == 0
                jacobfourier.stroke(j,i)=chy{i}(j,:)*jacobian.stroke(:,i);
            end
            jacobfourier.eqi(j,i)=chy{i}(j,:)*jacobian.eqi(:,i);
        end
    end

    if systemtype == 1

        if (isstruct(f))
            T = f.T;
            for fn = fieldnames(f).'
                p.(fn{:}) = f.(fn{:});
            end
        else
            p = makeGait(f);
            T = 2*pi/f(end,1);
        end

        s.costfunction = gopt.costfunction;

        % Get the gradient of cost based on inertia-dominated system
        jacobfourier.stroke = inertia_cost_gradient(s,n,f,T,p,'discrete');

        % With cost as time period to execute the gait, the gradient of
        % cost for inertial systems is the gradient of cost with period 1
        % divided by (4*T^3)
        jacobfourier.stroke = jacobfourier.stroke./(4*totalstroke^3);

        % Exclude the gradient of stroke for frequency terms,
        % since the inertia gradient includes those terms
        jacobfourier.stroke(end,:) = [];
    end
    % jacobianeqifourier commented out at Hossein's suggestion - don't
    % necessarily want the points to all be equally spaced
    jacobfourier.eff = jacobfourier.disp/totalstroke-...
        avgd_vel_dir*jacobfourier.stroke/totalstroke^2;%+jacobfourier.eqi;

end

