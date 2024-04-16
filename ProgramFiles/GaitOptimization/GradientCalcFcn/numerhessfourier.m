function hess = numerhessfourier(f,s,options)
    % By using central finite differenece, calculate the hessian of the
    % displacement and the cost. This function use the parallal for-loop,
    % and needs a parallal computing toolbox.
    %
    % Inputs:
    %   f: Fourier coefficients that parametrize the gait.
    % s: System file which contains the connection vector field, CCF's and
    %   metric data
    % dimension: Indicates the number of shape variables of the system
    % n: The number of points desired in a direct transcription parametrization
    %   of the gaits
    % direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
    % costfunction: costfunction type to optimize for

    dimension = options.dimension;

    dx = 1e-6;
    m = size(f,1)-1;

    % Initialize the struct to contain the intermediate jacobian.
    dgrad1 = cell(m,dimension);
    dgrad2 = cell(m,dimension);

    hess = struct();

    parfor i = 1:m
        for j = 1:dimension
            % give the perturbation to i,j-th entries of df.
            df=zeros(m+1,dimension);
            df(i,j) = dx;
            
            % calculate the jacobians.
            dgrad1{i,j} = evaluate_jacobian_fourier(f-df,s,options);
            dgrad2{i,j} = evaluate_jacobian_fourier(f+df,s,options);

        end
    end
    
    % calculate the hessians from the jacobians.
    hess.disp = cell(dimension,dimension);
    hess.stroke = cell(dimension,dimension);
    hess.eqi = cell(dimension,dimension);
    for i = 1:dimension
        for k = 1:m
            temphessdisp=(dgrad2{k,i}.disp - dgrad1{k,i}.disp).'/(2*dx);
            temphessstroke=(dgrad2{k,i}.stroke - dgrad1{k,i}.stroke).'/(2*dx);
            temphesseqi=(dgrad2{k,i}.eqi - dgrad1{k,i}.eqi).'/(2*dx);
            for j = 1:dimension
                hess.disp{i,j}(k,:) = temphessdisp(j,:);
                hess.stroke{i,j}(k,:) = temphessstroke(j,:);
                hess.eqi{i,j}(k,:) = temphesseqi(j,:);
            end
        end
    end
end

