function [A,Aeq,gA,gAeq]=nonlcon(y,s,gopt)
%%%%%%%%% 
%
%This function imposes the nonlinear constraint that all the points forming the gait stay within bounds
%
%Inputs:
%
%y: Fourier series coefficients that describe the gait
%s: System file which contains the connection vector field, CCF's and
%   metric data
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
%dimension: Indicates the number of shape variables of the system
%
%%%%%%%%%

n = gopt.npoints;
dimension = gopt.dimension;
if (gopt.issubopt)
    direction = gopt.subdirection;
else
    direction = gopt.direction;
end
constraint = gopt.constraint;
s.costfunction = gopt.costfunction;

%% Getting a discrete waypoint series of Fourier Gait.

% The first step is to obtain a direct transciption parametrization of 
% the gait from the fourier series parametrization
y1 = path_from_fourier(y,n-1,dimension);

% Since the shape boundary is the linear inequality,
% this function does not contain the inequality.
A = [];
gA = [];

%% Make sure the frequency doesn't get changed from 2*pi
Aeq = y(end,:)' - 2*pi;
if nargout > 2
    gAeq = zeros(numel(y),dimension);
end

%% The optimal pacing.
% [jacobfourier,~,~,eqi] = evaluate_jacobian_fourier(y,s,n,dimension,direction);
% Aeq(end+1) = eqi;
% if nargout > 2
%     geqi = jacobfourier.eqi;
%     geqi = reshape(geqi,[size(y,1)*dimension,1]);
%     gAeq(:,end+1) = geqi;
% end

%% Make sure that the displacement in the other two directions is zero
% Assign a time period for executing the gait
if any(constraint)
    [~, net_disp_opt,cost] = evaluate_displacement_and_cost1(s,y);
    net_disp_opt = se2log(net_disp_opt);
end

% Constrain solutions to only displace in the desired direction
if constraint(1) == 1
    for idx=1:3
        if idx ~= direction
            Aeq(end+1) = net_disp_opt(idx);
            if nargout > 2
                gopt_copy = gopt;
                gopt_copy.direction = idx;
                jacobfourier = evaluate_jacobian_fourier(y,s,gopt_copy);
                gdisp = [jacobfourier.disp; zeros(1,dimension)];              
                gdisp = reshape(gdisp,[size(y,1)*dimension,1]);
                gAeq(:,end+1) = gdisp.';
            end
        end
    end
end

% If optimizing for translation, restrict to zero rotation.
if constraint(2) == 1
    if direction ~=3
        Aeq(end+1) = net_disp_opt(3)/cost;
        if nargout > 2
            jacobfourier = evaluate_jacobian_fourier(y,s,n,dimension,3);
            geff = [jacobfourier.eff; zeros(1,dimension)];
            geff = reshape(geff,[size(y,1)*dimension,1]);
            gAeq(:,end+1) = geff.';
        end
    end
end

if length(constraint) >= 3
    if constraint(3) ~= 0
        Aeq(end+1) = net_disp_opt(direction) - constraint(3);
        if nargout > 2
            jacobfourier = evaluate_jacobian_fourier(y,s,n,dimension,direction);
            gdisp = [jacobfourier.disp; zeros(1,dimension)];
            gdisp = reshape(gdisp,[size(y,1)*dimension,1]);
            gAeq(:,end+1) = gdisp.';
        end
    end
end

% if(nargout > 3)
%     gradA = [];    
%     
%     gradAeq = ones(size(A,1),dimension);
%     
%     jacobiandisp = zeros(n,dimension);
%     for i=2:1:n-1
%         jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
%     end
%     jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
%     jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);
% end
% end
