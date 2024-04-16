function [jacobian,avgd_vel,totalstroke,avgd_vel_dir] = evaluate_jacobian(f,s,gopt)
%%%%%%%%%%%%%
% This function calculates the gradient of displacement and cost function
% with respect to the coefficients obtained by 
% the waypoint series parametrization
% 
% Inputs: 
% 
% f: Matrix containing the Fourier series coefficients
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
% jacobfourier: the struct which contains the gradient of the functions.
%   (jacobfourier.disp, jacobfourier.stroke, jacobfourier.eqi,
%   jacobfourier.total)
% lineint : the displacement for the gait
% totalstroke : the cost for the gait
%%%%%%%%%%%%%

n = gopt.npoints;
dimension = gopt.dimension;
if (gopt.issubopt)
    direction = gopt.subdirection;
else
    direction = gopt.direction;
end

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

%% Obtaining points from the gait parameters
if (isstruct(f))
    T = f.T;
    for fn = fieldnames(f).'
        p.(fn{:}) = f.(fn{:});
    end
    y = [f.phi_def{1}(linspace(0,T,n+1)).',f.phi_def{2}(linspace(0,T,n+1)).'];
    y = y(1:end-1,:);
else
    p = makeGait(f);
    y = path_from_fourier(f,n,dimension);
    T = 2*pi/f(end,1);
    y = y(1:end-1,:);
end

%% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
s.costfunction = gopt.costfunction;
[~, net_disp_opt, totalstroke] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');

% The coordinate of gradient is on the exponential coordinate, 
% but the net displacement derived by line integral is on the original
% coordinate. Take the matrix logarithms on the net displacement.
avgd_vel = se2log(net_disp_opt);

% displacement produced in the chosen direction produced on executing 
% the gait measured in the optimal coordinates 
avgd_vel_dir = avgd_vel(direction);

% If efficiency is negative, reversing the order of points so that
% efficiency is positive
% if lineint<0
%     lineint= -lineint;
%     y=flip(y);
%     invert=1;
% else
%     invert=0;
% end

%% Preliminaries for gradient calculation
[ccf,dccf] = ccf_interpolator(y,s,n,dimension,direction);
[metric,metricgrad] = metric_interpolator(y,s,n,dimension);

%% Jacobianstroke is the gradient of cost. 
% Contrigrad is the contribution to the gradient ue to the metric changing
jacobian = struct();

%% Checking if the system is inertia-dominated or drag-dominated
if systemtype == 0
    % Get the gradient of cost based on drag-dominated system
    [jacobian.stroke,l] = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
elseif systemtype == 1
    % The gradient of cost for inertia-dominated system should be
    % calculated at evaluate_jacobian_fourier.
    jacobian.stroke = [];
else
    error('Unexpected system type at cost gradient calculation!')
end

%% Jacobiandisp is the gradient of displacement.
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. It's input arguments are the coordinates of 
% the (i-1)th, ith and (i+1)th point, CCF value at point i and the dimension of     
% the shape space (dimension)
jacobian.disp = zeros(n,dimension);
for i=2:1:n-1
    jacobian.disp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension,dccf(i,:));
end
jacobian.disp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension,dccf(1,:));
jacobian.disp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension,dccf(n,:));

%% Jacobianeqi is the concentration gradient. 
%It is the term that keeps points eqi distant from each other and prevents crossover of gait.
if systemtype == 0
    jacobian.eqi = jacobianeqicalculator(y,n,dimension,metric);
else
    jacobian.eqi = zeros(n,dimension);
end

%% properly ordering gradients depending on whether lineint was negative
% if invert
%     y=flip(y);
%     jacobian.disp=flip(jacobian.disp);
%     jacobian.eqi=flip(jacobian.eqi);
%     jacobian.stroke=flip(jacobian.stroke);
% end

end