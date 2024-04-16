function yi=optimalgaitfamilygenerator(s,y0,mu0,stretch,handles,gopt)

    nfparam = gopt.nfparam;
    dimension = gopt.dimension;
    direction = gopt.direction;

    % Let user know starting the family optimizer.
    disp(" ");
    disp(" ");
    disp("------- Starting the family optimizer -------");
    disp("Finding the initial Lagrange multipliers.");

    if isfield(gopt,"subdirection")
        direction = [direction gopt.subdirection];
        anchor = gopt.anchor;
    end

    % Calculate the lagrange multipliers for initial configuration.
    [multp0,c0] = initlagrmultp(s,y0,mu0,gopt);

    % Convert the gait parameter to 1d array.
    y0 = reshape(y0(1:end-1,:),[],1);

    % Include lagrange multipliers and continuation variables in
    % optimization variables.
    y0 = [y0; multp0; c0];

    disp("Starting the numerical continuation methods.")

    if ~isfield(gopt,"subdirection")
        cont_fcn = @(y) stepcontfunc(y,s,gopt);
    else
        cont_fcn = @(y) steeringcontfunc(y,s,gopt);
    end

    ncp_fcn = @(y) noncomplfun(y,gopt);
    out_fcn = @(y,i) contoutfun(y,i,s,stretch,handles,gopt);
    options = contoptions('algorithm','PseudoArcLen','h',0.1,'nulltol',0.05,'ncpfcn',ncp_fcn,'outfcn',out_fcn);
    q = contmethod(cont_fcn,y0,options);

    yf = q(1:(nfparam-1)*dimension,:).';
   

%     choose the optimal gaits (1,0.75,0.5,0.25 of max-eff gait)
    s.costfunction = gopt.costfunction;
    net_disp_opt = zeros(3,size(yf,1));
    cost = zeros(1,size(yf,1));
    yi = cell(1,size(yf,1));
    for i = 1:size(yf,1)
        yi{i} = reshape_parameter_to_twod(yf(i,:).',gopt);
        [~,net_disp_opt(:,i),cost(i)] = evaluate_displacement_and_cost1(s,yi{i});
        net_disp_opt(:,i) = se2log(net_disp_opt(:,i));
    end

    if ~isfield(gopt,'subdirection')
        c = net_disp_opt(direction,:);
    else
        eff1 = net_disp_opt(direction(1),:)./cost;
        eff2 = net_disp_opt(direction(2),:)./cost;
        eff1 = (eff1 - anchor(2,1))/(anchor(1,1)-anchor(2,1));
        eff2 = (eff2 - anchor(1,2))/(anchor(2,2)-anchor(1,2));

        c = atan2(eff2,eff1);

        m = 25;
        newc = linspace(c(1),c(end),m);

        newq = (interp1(c.',q.',newc.','linear'));
        newq = [newq zeros(size(newq,1),1)];
        newq = num2cell(newq,2);

        gopt.neq = 2;
        nftotal = (nfparam-1)*dimension;

        options = contoptions('algorithm','PseudoArcLen','h',0.1,'nulltol',0.05,'ncpfcn',ncp_fcn,'outfcn',out_fcn);
        for i = 1:size(newq,1)
            [y0,~,~,mu0,~] = contvardistributor(newq{i}(1,:).',gopt);
            
            y0_2d = reshape_parameter_to_twod(y0,gopt);

            [multp0,c0,c02] = initlagrmultp(s,y0_2d,mu0,gopt);
            gopt.c02 = c02;
            y0 = [y0; multp0; c0];
            cont_fcn = @(y) steeringcontfunc(y,s,gopt);
            newq{i} = contmethod(cont_fcn,y0,options);
        end

        for i = 1:size(newq,1)
            for j = 1:size(newq{i},2)
                yi{i,j} = reshape_parameter_to_twod(newq{i}(1:nftotal,j),gopt);
            end
        end
    end

    ngaits = 4;

    cj = zeros(1,ngaits);
    [cj(1),j(1)] = max(c);
    j = zeros(1,ngaits);

    for i = 2:ngaits
        stepDisp = cj(1)*(ngaits+1-i)/ngaits;

        [~,j(i)] = min(abs(stepDisp-c));
        cj(i) = stepDisp - c(j(i));
    end

    y1=cell(ngaits,1);
    y=cell(ngaits,1);
    for i = 1:ngaits
        if i > 1
            if cj(i) > 0
                y11 = path_from_fourier(yi{j(i)-1},npoints,dimension);
                y12 = path_from_fourier(yi{j(i)},npoints,dimension);

                c1 = c(j(i)-1);
                c2 = c(j(i));
                cmid = c2 + cj(i);
            else
                y11 = path_from_fourier(yi{j(i)},npoints,dimension);
                y12 = path_from_fourier(yi{j(i)+1},npoints,dimension);
                c1 = c(j(i));
                c2 = c(j(i)+1);
                cmid = c1 + cj(i);
            end
            y1{i} = (y11*(c2-cmid) + y12*(cmid-c1))/(c2-c1);
        else
            y1{i} = path_from_fourier(yi{j(i)},npoints,dimension);
        end

        % path_from_fourier returns a self-connected gait, so remove the last point
        % to give what optimalgaitgenerator expects to return
        y1{i} = y1{i}(1:end-1,:);
        y{i}=y1{i}(:);
    end
end