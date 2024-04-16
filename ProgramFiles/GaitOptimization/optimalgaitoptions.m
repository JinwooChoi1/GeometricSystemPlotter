function options = optimalgaitoptions(s,options)
    arguments
        s struct
        options.dimension (1,1) double = 0
        options.npoints (1,1) double = 100
        options.nfourier (1,1) double = 4
        options.boundarydensity (1,1) double = 1
        options.direction (1,1) double
        options.subdirection (1,1) double
        options.issubopt (1,1) logical = false
        options.costfunction string = 'none'
        options.constraint
        options.familyoptimization logical = false
    end
    
    %%%%%%%%%%%%%%
    % A,b: There are lower and upper bound of shape variables for each point
    %   which is obtained from the grid inside which an optimal gait is
    %   desired. This relationship is expressed as A*p-b <= 0.
    %%%%%%%%%%%%%%

    %% Calculate nfparam. Total number of parameters per dimension.
    options.nfparam = (options.nfourier+1)*2;
    options.nftotal = options.nfparam*options.dimension;

    %% Set lower and upper bounds on optimization

    % Get the lower and upper bounds on the grid
    lvals = s.grid_range(1:2:end);
    uvals = s.grid_range(2:2:end);

    % Get the center of the grid
    mvals = (uvals + lvals)/2;

    % Get the lower and upper grid bounds relative to the center of the grid
    dlvals = lvals-mvals;
    duvals = uvals-mvals;

    % scale the differences by .95
    sdlvals = dlvals * 0.95;
    sduvals = duvals * 0.95;

    % Add the scaled differences to the center
    slvals = mvals+sdlvals;
    suvals = mvals+sduvals;

    % Tile these out to a matrix with as many rows as there are time points
    lb = 0.95 * repmat(slvals,[floor(options.npoints*options.boundarydensity),1]);
    ub = 0.95 * repmat(suvals,[floor(options.npoints*options.boundarydensity),1]);

    % Stack the lower bound values and the upper bound values
    options.lb = lb(:);
    options.ub = ub(:);

    
    % Set a linear inequality constraint for a lower and upper bound
    % A*y <= b
    dimension = options.dimension;

    nlbpoints = floor(size(options.lb,1)/dimension);
    
    y0 = zeros(options.nfparam,dimension);
    y0(end,:) = ones(1,dimension)*2*pi;

    chy = chy_generator(y0,nlbpoints,dimension);
    options.A = cell(dimension,dimension);
    for i = 1:dimension
        for j = 1:dimension
            if i == j
                options.A{i,j} = [-chy{i}.' zeros(nlbpoints,1)];
                options.A{i+dimension,j} = [chy{i}.' zeros(nlbpoints,1)];
            else
                options.A{i,j} = zeros(nlbpoints,options.nfparam);
                options.A{i+dimension,j} = options.A{i,j};
            end
        end
    end
    options.A = cell2mat(options.A);
    options.b = [-options.lb; options.ub];
    options.Aeq=[];
    options.beq=[];


    %% Set Objective and Constraint functions
    options = setgaitoptfcn(s,options);
end