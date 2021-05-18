function y=optimalgaitgenerator(s,dimension,npoints,a,lb,ub,stretch,direction,costfunction,handles)
%%%%%%%%%%%%%%
% This function takes an input gait and runs fmincon to find the neareast locally 
% optimal gait

%Inputs:
%
%s: System file which contains the connection vector field, CCF's and
%   metric data
%dimension: Indicates the number of shape variables of the system
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
% a: n x m array, for n waypoints in m dimensions
% lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% direction: Direction to optimize travel: 1=x,2=y,3=theta
% costfunction: costfunction type to optimize for
%           Options: Path length (metric), torque squared (metric),
%           covariant acceleration (metric), path length (coordinate),
%           acceleration (coordinate)
% 
% 
% Outputs: 
% 
% y: Matrix whose values indicate coordinates of points which form the optimal gait
%%%%%%%%%%%%

% n=npoints;
% for i=1:1:npoints
%     a1(1,i)=1*cos(2*pi*(i-1)/n);
%     a2(1,i)=1*cos(2*pi*(i-1)/n+pi/2);
% end

% n=npoints;
% P1(:,1)=a1(1,1:n)';
% P1(:,2)=a2(1,1:n)';
% P1(end+1,:) = P1(1,:); % Close the loop on the gait

% For minimal refactoring, mapping a -> P1
P1 = a;

% % Close the loop of the gait if necessary
% if P1(end,:) ~= P1(1,:)
%     P1(end+1,:) = P1(1,:);
% end


%% Finding fourier coeffecients.
% The first step is to go from a direct transcription of the initial gait
% to a fourier based parametrization. 
% fa is a cell where the ith element contains the coefficients for the fourier parametrization along the ith direction 

% Time period of gait is 1 second in order to handle calculations performed
% for inertial gait optimization
t = linspace(0,1,size(P1,1));

fa=cell(dimension);
% The bounds ensure the fourier series terms have the right period
 % Bound only the frequency to be 2*pi, such that period = 1
options = fitoptions('fourier4');
options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi];
options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi];

for i=1:1:dimension
    fa{i}=fit(t',P1(:,i),'fourier4',options);
end

%% The next step is to setup the fmincon call. 
% y0 is the marix of all fourier series coefficients that describe the
%   initial gait
% nonlcon is the function that imposes the constraints that all the points
%   stay inside the grid
% outfun is the function that updates the gait on the GUI after every iteration 

A=[];
b=[];
Aeq=[];
beq=[];

nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};%
 lb1=[];
 ub1=[];

 y0 = zeros(length(nu),dimension);
for i=1:dimension
    for j=1:length(nu)
        y0(j,i)=fa{i}.(nu{j});
    end
end

writerObj = [];
% % Uncomment this section if you'd like to record a video of the
% % optimizer's steps in the shape space
% writerObj = VideoWriter('cost_as_time_period_circle_start.mp4','MPEG-4');
% writerObj.FrameRate = 5;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% figure(5);
% subplot(1,2,1)
% contour(s.grid.eval{1},s.grid.eval{2},s.DA_optimized{1},'LineWidth',1.5)
% axis square
% hold on

s.costfunction = costfunction;
global bestCost bestDisp bestEff;
bestEff = 0;

%Suppress warning for annoying thing in jacobianeqicalculator
warning('off','MATLAB:sqrtm:SingularMatrix');

try
 options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display','iter','Algorithm','sqp','CheckGradients',false,'FiniteDifferenceType','central','MaxIter',4000,'MaxFunEvals',20000,'TolCon',10^-3,'StepTolerance',1e-6,'OutputFcn', @(y,optimValues,state) outfun(y,optimValues,state,stretch,s,handles));
catch
   error('This code requires the global optimization toolbox to run') 
end
 [yf, ~,~,~]=fmincon(@(y) solvedifffmincon(y,s,npoints,dimension,direction,lb,ub,writerObj),y0,A,b,Aeq,beq,lb1,ub1,@(y) nonlcon(y,s,npoints,dimension,lb,ub),options);

% % Uncomment this if you uncommented the section above so that the video
% % writer object is closed appropriately.
% close(writerObj);

printstuff = 1;
if printstuff
    disp(['Optimal Efficiency: ',num2str(bestEff)]);
    disp(['Optimal Displacement: ',num2str(bestDisp)]);
    disp(['Optimal Cost: ',num2str(bestCost)]);
end

%% Getting point position values from the result of fmincon
% This section helps us go back to a direct transcription parametrization
% of the optimal gait from a fourier series parametrization. y is a column vector
% that contains coordinates of all points forming the optimized gait

y1 = path_from_fourier(yf,npoints,dimension);
% path_from_fourier returns a self-connected gait, so remove the last point
% to give what optimalgaitgenerator expects to return
y1 = y1(1:end-1,:);
y=y1(:);

%% Uncomment for plotting the optimized gait. Potentially useful while debugging.
% for i=1:n
%     xf(i)=y(i);
%     yf(i)=y(n+i);
% end
% % 
% % for i=1:n`
% %     xf(i)=P1(i,1);
% %     yf(i)=P1(i,2);
% % end
% % 
% figure(11)
% hold on
% plot(xf,yf,'-o')

end

function [f,g]=solvedifffmincon(y,s,n,dimension,direction,~,~,~)%,lb,ub,writerObj)
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
% lb: Lower bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside 
%   which an optimal gait is desired
% 
% Outputs:
% 
% f: Objective function value (This is negative of efficiency by default, can be
%   changed to displacement)
% g: Gradient of the objective function
%%%%%%%%%%%%%

%% Obtaining points from fourier coefficients
% The first step is to obtain a direct transcription of the gait from the
% fourier series parametrization. y1 is the matrix of coordinate values of
% the points forming the gait.
afactor=0.001;
coeff=y;
y1 = path_from_fourier(y,n,dimension);
y1 = y1(1:end-1,:); % Remove the last points because path_from_fourier returns self-connected gait

%% Calculating cost and displacement per gait

w1 = y(end,1); % Frequency of Fourier transform
w2 = y(end,2);
% Assign a time period for executing the gait
T = 2*pi/w1;

% Define phi_def = [alpha1, alpha2] as a function of time t such that the
% array returns the shape variables given by the fourier coefficients at a
% time t
p = makeGait(y);
        
% % Uncomment this section to verify that the shape variables and derivatives
% % have been derived appropriately
% valid = verify_shape_equations(p);
% valid_M = verify_mass_derivative(s);
% validate_cost_gradient(s,n,y,T,p);

% Reassignment of point transcription to variable y for compatibility with
% old version
clear y
y=y1;
clear y1;

global bestCost bestDisp bestEff;

% Calculate displacement, cost and efficiency of a gait
% Note that, for inertial cost, cost is returned as the integral of torque
% squared, while for drag-based systems, cost is the path length of the
% gait
[~, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');
lineint=net_disp_opt(direction); % displacement produced in the chosen direction produced on executing the gait measured in the optimal coordinates 

% Assign value for totalstroke, i.e. the cost metric used for efficiency
% calculation
if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
    % With cost as time period, period of the gait is the cost to execute
    % the gait at unit torque squared to the 1/4th power
    totalstroke = cost^(1/4);
else
    % Drag systems have cost as path length, so no modification is needed
    totalstroke = cost;
end

% If efficiency is negative, reversing the order of points so that
% efficiency is positive
if lineint<0
    lineint= -lineint;
    ytemp=y;
    for i=1:n
        y(i,:)=ytemp(n+1-i,:);
    end
    invert=1;
else
    invert=0;
end

%% Preliminaries for gradient calculation
% Preallocating memory for variables which we will need in further
% calculation 
yvalues=cell(n,dimension); % Cell representation of the coordinates of all points forming the gait
interpstateccf=cell(1,dimension); % Variable which will store the ccf function grid used for interpolation
interpmetricgrid=cell(1,dimension);  % Variable which will store the metric grid used for interpolation
ccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store ccf function at each point
metric1=zeros(n,dimension,dimension); % Variable which will store metric at each point in the form of a matrix
metric = repmat({zeros(dimension)},[n 1]); % Variable which stores the metric at each point in the form of a 2x2 matrix
metricgrad1=zeros(n,dimension,dimension,dimension); % Variable which will store gradient of metric at each point in the form of a matrix
metricgrad = repmat({zeros(dimension)},[n,dimension]); % Variable which will store gradient of metric at each point in the form of a matrix

% Interpolation to calculate all the variables needed for gradient
% calculation
for i=1:1:n
    for j=1:1:dimension
        yvalues{i,j}=y(i,j);
    end
end

y_for_interp = mat2cell(y,size(y,1),ones(1,size(y,2)));

for j=1:1:dimension
    interpstateccf{j}=s.grid.eval{j,1};
    interpmetricgrid{j}=s.grid.metric_eval{j,1};
end

for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y_for_interp{:},'spline');
end

for j=1:1:dimension
    for k=1:1:dimension
        metric1(:,j,k)=interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y_for_interp{:},'spline');
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:n
       metric{i}=eye(dimension);
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    for i=1:n
        for j=1:1:dimension
           for k=1:1:dimension
               metric{i}(j,k)=metric1(i,j,k);
           end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for i = 1:n
        metric{i} = metric{i}*metric{i};
    end
end

if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'acceleration coord')
    for i=1:dimension
        for j=1:n
            metricgrad{j,i} = zeros(dimension);
        end
    end
elseif strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
    y2 = zeros(size(y));
    y1 = y2;
    for l=1:1:dimension
        for m=1:1:dimension
            if m==l
               y2(:,m)=y(:,m)+afactor*ones(length(y),1);
               y1(:,m)=y(:,m)-afactor*ones(length(y),1);
            else
               y2(:,m)=y(:,m);
               y1(:,m)=y(:,m);
            end
        end
        y2_for_interp = mat2cell(y2,size(y,1),ones(1,size(y,2)));
        y1_for_interp = mat2cell(y1,size(y,1),ones(1,size(y,2)));
        for j=1:1:dimension
            for k=1:1:dimension
                metricgrad1(:,l,j,k)=(interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y2_for_interp{:},'spline')...
                    -interpn(interpmetricgrid{:},s.metricfield.metric_eval.content.metric{j,k},y1_for_interp{:},'spline'))/(2*afactor);
            end
        end
        for i=1:n
            for j=1:1:dimension
                for k=1:1:dimension
                    metricgrad{i,l}(j,k)=metricgrad1(i,l,j,k);
                end
            end
        end
    end
end
if strcmpi(s.costfunction,'pathlength metric2')
    for l = 1:dimension
        for i = 1:n
            metricgrad{i,l} = metricgrad{i,l}*metric{i}+metric{i}*metricgrad{i,l};
        end
    end
end

%% changey/dcoeff tells us how much each point moves when a fourier series variable is changed
% chy is a cell with as many entries as the dimension of the shape space
% ith element of chy is a matrix where the (j,k)th entry tells us the change in the ith coordinate
% of the kth point of the gait resulting from a unit change in the jth
% fourier coefficient corresponding to the ith dimension of the shape space

chy=cell(dimension,1);
% Create vector of time values at which to evaluate points of gait
t = linspace(0,T,n);
for i=1:1:dimension
    for j=1:1:n
        chy{i}(:,j)=[1;cos(t(j)*coeff(end,i));sin(t(j)*coeff(end,i));cos(2*t(j)*coeff(end,i));sin(2*t(j)*coeff(end,i));cos(3*t(j)*coeff(end,i));sin(3*t(j)*coeff(end,i));cos(4*t(j)*coeff(end,i));sin(4*t(j)*coeff(end,i))];%cos(5*t(j)*coeff(end,i));sin(5*t(j)*coeff(end,i))];%;cos(6*t(j)*coeff(end,i));sin(6*t(j)*coeff(end,i))];%
    end
end

%% Jacobianstroke is the gradient of cost. 
%Contrigrad is the contribution to the gradient due to the metric changing
switch s.costfunction
    case {'pathlength metric','pathlength coord','pathlength metric2'}
        % Get the gradient of cost based on drag-dominated system
        jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad);
    case {'torque','covariant acceleration','acceleration coord','power quality'}
        % Get the gradient of cost based on inertia-dominated system
        inertia_cost_grad = inertia_cost_gradient(s,n,coeff,T,p,'discrete');
        
        % With cost as time period to execute the gait, the gradient of
        % cost for inertial systems is the gradient of cost with period 1
        % divided by (4*T^3)
        inertia_cost_grad = inertia_cost_grad./(4*totalstroke^3);
    otherwise
        error('Unexpected system type at cost gradient calculation!')
end

%% Jacobiandisp is the gradient of displacement.
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. It's input arguments are the coordinates of 
% the (i-1)th, ith and (i+1)th point, CCF value at point i and the dimension of     
% the shape space (dimension)

jacobiandisp = zeros(n,dimension);
for i=2:1:n-1
    jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
end
jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);

%% Jacobianeqi is the concentration gradient. 
%It is the term that keeps points eqi distant from each other and prevents crossover of gait.

jacobianeqi = jacobianeqicalculator(y,n,dimension,metric);

%% properly ordering gradients depending on wether lineint was negative or positive
if invert
        jacobiandisptemp=jacobiandisp;
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
                jacobianstroketemp=jacobianstroke;
        end
        jacobianeqitemp=jacobianeqi;
    for i=1:1:n
        jacobiandisp(i,:)=jacobiandisptemp(n+1-i,:);
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
            jacobianstroke(i,:)=jacobianstroketemp(n+1-i,:);
        end
        jacobianeqi(i,:)=jacobianeqitemp(n+1-i,:);
    
    end
end

%% fourier series version of all gradients

% We then obtain gradients in a fourier series parametrization by
% projecting the gradients from the direct transcription space onto the
% fourier coefficient space
for i=1:1:dimension
    for j=1:1:9 
        jacobiandispfourier(j,i)=chy{i}(j,:)*jacobiandisp(:,i);
        if strcmpi(s.costfunction,'pathlength coord') || strcmpi(s.costfunction,'pathlength metric') || strcmpi(s.costfunction,'pathlength metric2')
            jacobianstrokefourier(j,i)=chy{i}(j,:)*jacobianstroke(:,i);
        end
        jacobianeqifourier(j,i)=chy{i}(j,:)*jacobianeqi(:,i);
    end
end
if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
    % Inertia cost gradient is already in terms of the fourier coefficients
    jacobianstrokefourier = inertia_cost_grad;
    % Add zero terms to the gradient of displacement and jacobianeqi for
    % frequency terms, since the inertia gradient includes those terms
    jacobiandispfourier = [jacobiandispfourier;zeros(1,dimension)];
    jacobianeqifourier = [jacobianeqifourier;zeros(1,dimension)];
    % Calculate the total gradient of efficiency
    % jacobianeqifourier commented out at Hossein's suggestion - don't
    % necessarily want the points to all be equally spaced
    totaljacobianfourier = jacobiandispfourier/totalstroke-lineint*jacobianstrokefourier/totalstroke^2;%+jacobianeqifourier;
    % reset the gradient of the frequency terms to be zero so they aren't
    % changed
    totaljacobianfourier(end,:) = zeros(1,dimension);
else
    totaljacobianfourier = jacobiandispfourier/totalstroke-lineint*jacobianstrokefourier/totalstroke^2+jacobianeqifourier;
end

%% minimizing negative of efficiency(or displacement)
f=-lineint/(totalstroke); % Optimizing for displacement over cost
%     f = -lineint; % Optimizing for displacement only

if abs(f) > bestEff
    bestEff = abs(f);
    bestDisp = abs(lineint);
    bestCost = abs(totalstroke);
end

if nargout>1
    if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
        % Return the gradient of efficiency as previously calculated for
        % inertia systems
        g = -totaljacobianfourier; % Optimizing for displacement over cost
%         g = -jacobiandispfourier; % Optimizing for displacement only
    else
        % Return the gradient of efficiency plus row of zeros for frequency
        % terms for drag systems
        g=[-totaljacobianfourier;zeros(1,dimension)]; % Optimizing for displacement over cost
%         g = [-jacobiandispfourier;zeros(1,dimension)]; % Optimizing for displacement only
    end
end

end

%Makes gait definition structure of cell functions from
%fourier coefficients y for arbitrary dimension.
function gait = makeGait(y)

    ndim = size(y,2);
    
    gait.phi_def = cell(1,ndim);
    gait.dphi_def = cell(1,ndim);
    gait.ddphi_def = cell(1,ndim);
    
    for dim = 1:ndim
        
        w = y(end,dim);
        gait.phi_def{dim} = @(t) y(1,dim)+y(2,dim)*cos(w*t)+y(3,dim)*sin(w*t)+y(4,dim)*cos(2*w*t)+...
                                +y(5,dim)*sin(2*w*t)+y(6,dim)*cos(3*w*t)+y(7,dim)*sin(3*w*t)+...
                                +y(8,dim)*cos(4*w*t)+y(9,dim)*sin(4*w*t);
        gait.dphi_def{dim} = @(t) -w*y(2,dim)*sin(w*t)+w*y(3,dim)*cos(w*t)-2*w*y(4,dim)*sin(2*w*t)+...
                                  +2*w*y(5,dim)*cos(2*w*t)-3*w*y(6,dim)*sin(3*w*t)+3*w*y(7,dim)*cos(3*w*t)+...
                                  -4*w*y(8,dim)*sin(4*w*t)+4*w*y(9,dim)*cos(4*w*t);
        gait.ddphi_def{dim} = @(t) -w^2*y(2,dim)*cos(w*t)-w^2*y(3,dim)*sin(w*t)-4*w^2*y(4,dim)*cos(2*w*t)+...
                                   -4*w^2*y(5,dim)*sin(2*w*t)-9*w^2*y(6,dim)*cos(3*w*t)-9*w^2*y(7,dim)*sin(3*w*t)+...
                                   -16*w^2*y(8,dim)*cos(4*w*t)-16*w^2*y(9,dim)*sin(4*w*t);
        
    end
                
end

%Finds value of gait cell function at given time
function state = readGait(gaitFun,t)

    ndim = numel(gaitFun);
    state = zeros(1,ndim);
    
    for i = 1:ndim
        state(i) = gaitFun{i}(t);
    end

end

function jacobianstroke = jacobianstrokecalculator(y,n,dimension,metric,metricgrad)
    % Calculates the gradient of cost for drag dominated systems.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    %   metricgrad: Gradient of Riemannian metric
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    delp = cell(1,n);
    for i=1:(numel(delp)-1)
        delp{i}=y(i+1,:)-y(i,:); % delp{i} is the vector joining the (i+1)th point to the ith point 
    end
    delp{end}=y(1,:)-y(n,:);

    jacobianstroke = zeros(n,dimension);
    contrigrad=zeros(n,dimension);
    for i=2:n-1
        for j=1:dimension
            %Contrigrad is the contribution to the gradient due to the metric changing
            contrigrad(i,j)=0.5*delp{i}*metricgrad{i,j}*delp{i}'/(2*l(i))+0.5*delp{i-1}*metricgrad{i,j}*delp{i-1}'/(2*l(i-1)); 
        end
        % Total gradient is the result of distance changing due to movement of point and the metric changing due to movement of the point
        jacobianstroke(i,:)=(-(((metric{i}+metric{i+1})/2)*delp{i}')'-(delp{i}*((metric{i}+metric{i+1})/2)))/(2*l(i))+...
            +((((metric{i-1}+metric{i})/2)*delp{i-1}')'+(delp{i-1}*((metric{i}+metric{i-1})/2)))/(2*l(i-1))+contrigrad(i,:); 
    end

    % Calculation for the 1st point and last point have to be done outside the
    % loop as the (i+1)th point for the last point is the first point and
    % (i-1)th point for the first point is the last point
    for j=1:dimension
        contrigrad(1,j)=0.5*delp{1}*metricgrad{1,j}*delp{1}'/(2*l(1))+0.5*delp{n}*metricgrad{1,j}*delp{n}'/(2*l(n));
    end
    jacobianstroke(1,:)=(-(((metric{1}+metric{2})/2)*delp{1}')'-(delp{1}*((metric{1}+metric{2})/2)))/(2*l(1))+...
        +((((metric{n}+metric{1})/2)*delp{n}')'+(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+contrigrad(1,:);

    for j=1:dimension
        contrigrad(n,j)=0.5*delp{n}*metricgrad{n,j}*delp{n}'/(2*l(n))+0.5*delp{n-1}*metricgrad{n,j}*delp{n-1}'/(2*l(n-1));
    end
    jacobianstroke(n,:)=(-(((metric{n}+metric{1})/2)*delp{n}')'-(delp{n}*((metric{n}+metric{1})/2)))/(2*l(n))+...
        +((((metric{n}+metric{n-1})/2)*delp{n-1}')'+(delp{n-1}*((metric{n}+metric{n-1})/2)))/(2*l(n-1))+contrigrad(n,:);
end

function jacobianeqi = jacobianeqicalculator(y,n,dimension,metric)
    % Calculates the gradient of the force driving the points along the
    % gait to be equally spaced.
    % Inputs:
    %   y: matrix containing points that transcribe the gait
    %   n: number of points in gait transcription
    %   dimension: number of shape variables
    %   metric: Riemannian metric
    
    jacobianeqi = zeros(n,dimension);
    
    %l is the vector containing metric weighted distances between neighbouring
    %points
    l = zeros(1,n);
    for i=1:(numel(l)-1)
        l(i)=sqrt((y(i+1,:)-y(i,:))*((metric{i}+metric{i+1})/2)*(y(i+1,:)-y(i,:))');
    end
    l(end)=sqrt((y(1,:)-y(n,:))*((metric{n}+metric{1})/2)*(y(1,:)-y(n,:))');
    
    for i=2:n-1
        len=sqrt((y(i+1,:)-y(i-1,:))*((metric{i-1}+metric{i+1})/2)*(y(i+1,:)-y(i-1,:))'); % metric weighted length between point (i-1) and (i+1)
        midpoint=y(i-1,:)+((y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2))/2; % location of midpoint of the line segment joining point (i-1) and (i+1)
        betacos=(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*((y(i,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i})/2))'/(l(i-1)*len);
        xhat=y(i-1,:)+(y(i+1,:)-y(i-1,:))*sqrtm((metric{i-1}+metric{i+1})/2)*l(i-1)*betacos/len; %projection of ith point onto the line joining the (i-1)th and (i+1)th points
        jacobianeqi(i,:)=midpoint-xhat; % gradient of the ith point is equal to the difference between the midpoint and the projection of ith point
    end

    len=sqrt((y(2,:)-y(n,:))*((metric{2}+metric{n})/2)*(y(2,:)-y(n,:))');
    midpoint=y(n,:)+((y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2))/2;
    betacos=(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*((y(1,:)-y(n,:))*sqrtm((metric{n}+metric{1})/2))'/(l(n)*len);
    xhat=y(n,:)+(y(2,:)-y(n,:))*sqrtm((metric{n}+metric{2})/2)*l(n)*betacos/len;
    jacobianeqi(1,:)=midpoint-xhat;

    len=sqrt((y(1,:)-y(n-1,:))*((metric{1}+metric{n-1})/2)*(y(1,:)-y(n-1,:))');
    midpoint=y(n-1,:)+((y(1,:)-y(n-1,:))*sqrtm((metric{1}+metric{n-1})/2))/2;
    betacos=(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*((y(n,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{n})/2))'/(l(n-1)*len);
    xhat=y(n-1,:)+(y(1,:)-y(n-1,:))*sqrtm((metric{n-1}+metric{1})/2)*l(n-1)*betacos/len;
    jacobianeqi(n,:)=midpoint-xhat;
end

function a=jacobiandispcalculator3(p1,p2,p3,ccf,dimension)
%%%%%%%%%
%
% jacobiandispcalculator3 is the function that calculates the gradient of 
% displacement for the ith point. 
% Its input arguments are the coordinates of the (i-1)th, ith and (i+1)th point,     
% CCF value at point i(ccf) and the dimension of the shape space (dimension)
%
%%%%%%%%%

l1=0; % variable for calculating length of the line segment joining the (i-1)th point with the (i+1)th point
base = zeros(1,dimension);
for i=1:numel(base)
    l1=l1+(p1(i)-p3(i))^2;
    base(1,i)=p3(i)-p1(i); % vector connecting the (i-1)th point and (i+1)th point  
end
%l=sqrt(l1); % length of the line segment joining the (i-1)th point with the (i+1)th point

jacobian = zeros(1,dimension);
for i=1:dimension
%    jacobian(1,i)=0;
    perp1=zeros(1,dimension);
    perp1(i)=1;
    %parcomp=base*perp1'/norm(base);
    %perp1-parcomp*base/norm(base);  %%recheck again
    perp=perp1;% Unit vector along the ith direction
    % The for loop below calculates the gradient along the ith direction by
    % treating the CCF as 2 forms. A specific (j,k) represents a component of the 2 form 
    for j=1:dimension-1 
        for k=1:dimension-j
            veca=zeros(1,dimension);
            vecb=zeros(1,dimension);
            veca(j)=1;
            vecb(j+k)=1;
            f=(j-1)*dimension-(j*(j-1))/2+k;
            jacobian(1,i)=jacobian(1,i)+0.5*ccf(f)*((veca*perp')*(vecb*base')-(vecb*perp')*(veca*base'));
        end
    end
end

a=jacobian;

end

function [A,Aeq]=nonlcon(y,s,n,dimension,lb,ub)
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
%lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
%ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% 
%%%%%%%%%

% % The first step is to obtain a direct transciption parametrization of the gait from the 
% % fourier series parametrization
y1 = path_from_fourier(y,n,dimension);
y2=y1(:);

%b=length(y2);

% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A = [A1;A2];


% Make sure the frequency doesn't get changed from 2*pi
Aeq = y(end,:) - 2*pi;

end


function y = path_from_fourier(f,n,dimension)
% Returns the shape space parametrization of the gait at n points when provided with
% the fourier coefficients f. The gait is returned as a self-closed gait
% (i.e. the first and last rows of the output y are the same point).
% Inputs:
%   f: Fourier coefficients that parametrize the gait.
%   n: Number of points that should compose the gait, less the final
%       self-closed point
%   dimension: Number of shape variables for system

    y = zeros(n+1,dimension);
    % Determine time period based on value of fourier frequency
    w = f(end,1);
    T = 2*pi/w;
    % Create time vector at which to evaluate points of gait
    t = linspace(0,T,n+1);
    % Evaluate the shape-space parametrization of the gait at every time
    % value in t
    for j=1:dimension
        for i=1:1:n+1
            y(i,j)=f(1,j)+f(2,j)*cos(w*t(i))+f(3,j)*sin(w*t(i))+f(4,j)*cos(2*w*t(i))+...
                +f(5,j)*sin(2*w*t(i))+f(6,j)*cos(3*w*t(i))+f(7,j)*sin(3*w*t(i))+...
                +f(8,j)*cos(4*w*t(i))+f(9,j)*sin(4*w*t(i));
        end
    end
end

function stop=outfun(y,optimValues,state,stretch,s,handles)
%%%%%%%%% 
%
%This function plots the current state of the gait on the sysplotter GUI
%after every iteration of the optimizer
%
%%%%%%%%% 

n=100;
dimension=length(y(1,:));

% % The if else statement below deletes gaits 2 iterations after they have been plotted
% if optimValues.iteration>2
%     children=get(gca,'children');
%     delete(children(6:10));
% else
% end

for thisAxes = [1:numel(handles.plot_thumbnails.Children)]
    
    axes(handles.plot_thumbnails.Children(thisAxes));

    % The if else statement below fades the gait plotted during the previous iteration
    if optimValues.iteration>1
        children=get(gca,'children');
        for idx = 1:numel(children)

            if iscell(children(idx).UserData) && strcmp(children(idx).UserData{1},'OptimizeTracer')
                children(idx).UserData = {'OptimizeTracer', children(idx).UserData{2}-1};

                if children(idx).UserData{2} == 0

                    delete(children(idx));

                else

                    children(idx).Color=[0.5 0.5 0.5];
                    children(idx).LineWidth=4;
                end
            end

        end
    %     children(1).Color=[0.5 0.5 0.5];
    %     children(2).Color=[0.5 0.5 0.5];
    %     children(3).Color=[0.5 0.5 0.5];
    %     children(4).Color=[0.5 0.5 0.5];
    %     children(5).Color=[0.5 0.5 0.5];
    % 
    %     children(1).LineWidth=4;
    %     children(2).LineWidth=4;
    %     children(3).LineWidth=4;
    %     children(4).LineWidth=4;
    %     children(5).LineWidth=4;
    else
    end

    % The if else statement below plots the gait after every iteration
    if optimValues.iteration>0
        y1 = path_from_fourier(y,n,dimension);
        hold on
        if stretch
            stretchnames = {'stretch','surface'};
            stretchname = stretchnames{stretch};

            [x_temp,y_temp,z_temp] = s.convert.(stretchname).old_to_new_points(y1(:,1),y1(:,2));
        else
            x_temp = y1(:,1);
            y_temp = y1(:,2);
            z_temp = zeros(size(y1(:,1)));
        end
        handle1=line('XData',x_temp,'YData',y_temp,'ZData',z_temp,'color','k','linewidth',3,'UserData',{'OptimizeTracer',2}); %#ok<NASGU>
        %plot_dir_arrows(y1(:,1),y1(:,2),2,'Color',[0 0 0],'LineWidth',3);
    else
    end

end

pause(0.05)
stop=false;
end


% Evaluate the body velocity and cost velocity (according to system metric)
% at a given time
function [gcirc, dcost] = get_velocities(t,s,gait,ConnectionEval,A,metric,dM)

	% Get the shape and shape derivative at the current time
    shape = zeros(size(s.grid.eval));
    dshape = zeros(size(s.grid.eval));
    ddshape = zeros(size(s.grid.eval));
    
	shape_gait_def = readGait(gait.phi_def,t);
	dshape_gait_def = readGait(gait.dphi_def,t);
    ddshape_gait_def = readGait(gait.ddphi_def,t);
    
    actual_size = min(numel(shape),numel(shape_gait_def));
    shape(1:actual_size) = shape_gait_def(1:actual_size);
    dshape(1:actual_size) = dshape_gait_def(1:actual_size);
    ddshape(1:actual_size) = ddshape_gait_def(1:actual_size);
  
    M_a = metric;
    
	shapelist = num2cell(shape);
	
    % If doing functional eval of system (not recommended)
	% Get the local connection and metric at the current time, in the new coordinates
	if strcmpi(ConnectionEval,'functional')
			
        A = s.A_num(shapelist{:})./s.A_den(shapelist{:});

        switch s.system_type
            case 'drag'
                metric = s.metric(shapelist{:});
            case 'inertia'
                error('Functional ConnectionEval method not supported for inertia systems!')
        end

    end
	
	% Get the body velocity at the current time
	%t;
    gcirc = - A * dshape(:);

    switch s.costfunction
        case {'pathlength metric','pathlength coord'}
            dcost = sqrt(dshape(:)'*metric*dshape(:));
        case 'pathlength metric2'
            dcost = sqrt(dshape(:)'*metric*metric*dshape(:));
        case 'torque'
            dcost = torque_cost(M_a,dM,shape,dshape,ddshape,metric);
        case 'covariant acceleration'
            dcost = acceleration_cost(M_a,dM,shape,dshape,ddshape,metric);
        case 'acceleration coord'
            dcost = ddshape(:)'*metric*ddshape(:);
        case 'power quality'
            dcost = power_quality_cost(M_a,dM,shape,dshape,ddshape);
    end
	
end

function dcost = torque_cost(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is torque squared.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM,shape,dshape);
    % Calculate the torque for this instant of time and return the inner
    % product of the torque with itself
    dtau = M*ddshape(:) + C;
    dcost = dtau'*dtau;
end

function dcost = acceleration_cost(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is covariant acceleration.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system
    C = calc_coriolis_matrix(dM,shape,dshape);
    cov_acc = ddshape(:) + inv(M)*C;
    dcost = cov_acc'*metric*cov_acc;

end

function dcost = power_quality_cost(M,dM,shape,dshape,ddshape)
% Calculates the incremental cost for an inertial system where cost is power quality.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM,shape,dshape);
    % Calculate the torque for this instant of time 
    dtau = M*ddshape(:) + C;
    % Calculate power quality
    dcost = (dshape(:)'*dtau)^2 - ((dshape(:)').^2*dtau.^2);
    dcost = dcost + 100;
end

function [net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,tspan,ConnectionEval,IntegrationMethod,resolution)
% Evaluate the displacement and cost for the gait specified in the
% structure GAIT when executed by the system described in the structure
% S.
%
% S should be a sysplotter 's' structure loaded from a file
% sysf_FILENAME_calc.mat (note the _calc suffix)
%
% P should be a structure with fields "phi_def" and "dphi_def", returning a
% vector of shapes and shape velocities respectively. If it is not
% convenient to analytically describe the shape velocity function,
% gait.dphi should be defined as 
%
% p.dphi =  @(t) jacobianest(@(T) p.phi (T),t)
%
% as is done automatically by sysplotter, but note that this will be slower
% than specifying a function directly
%
% ConnectionEval can specify whether the local connection should be generated from
% its original function definiton, or by interpolation into the evaluated
% matrix, but is optional. Valid options are 'functional' or
% 'interpolated'. Defaults to "interpolated", which significantly faster
% when calculating the local connection or metric from scratch takes
% appreciable computational time
%
% IntegrationMethod can specify whether ODE45 or a fixed point
% (euler-exponential) integration method should be employed. Defaults to
% fixed point, to reduce interpolation overhead computational times.
%
% RESOLUTION specifies the number of points for fixed-point resolution
% evaluation. A future option may support autoconvergence, but ODE
% performance with interpolated evaluation appears to be fast enough that
% refinement of fixed-point performance is on hold.
	

	% if no ConnectionEval method is specified, default to interpolated
	if ~exist('ConnectionEval','var')
		ConnectionEval = 'interpolated';
	end
    
    % if no IntegrationMethod is specified, default to ODE
	if ~exist('IntegrationMethod','var')
		IntegrationMethod = 'fixed_step';
	end

    % if no resolution is specified, default to 100 (this only affects
    % fixed_step integration)
	if ~exist('resolution','var')
		resolution = 100;
	end

    
    
	switch IntegrationMethod
		
		case 'fixed_step'
			
			[net_disp_orig, cost] = fixed_step_integrator(s,p,tspan,ConnectionEval,resolution);
        
        case 'ODE'

            % Calculate the system motion over the gait
            sol = ode45(@(t,y) helper_function(t,y,s,p,ConnectionEval),tspan,[0 0 0 0]');

            % Extract the final motion
            disp_and_cost = deval(sol,tspan(end));

            net_disp_orig = disp_and_cost(1:3);
            cost = disp_and_cost(end);
            
        otherwise
			error('Unknown method for integrating motion');
	end

	
	% Convert the final motion into its representation in optimal
	% coordinates
    startshape = zeros(size(s.grid.eval));
	startshape_def = readGait(p.phi_def,0);
        
    actual_size = min(numel(startshape),numel(startshape_def));
    startshape(1:actual_size) = startshape_def(1:actual_size);
    
	startshapelist = num2cell(startshape);
	beta_theta = interpn(s.grid.eval{:},s.B_optimized.eval.Beta{3},startshapelist{:},'spline');
	net_disp_opt = [cos(beta_theta) sin(beta_theta) 0;...
		-sin(beta_theta) cos(beta_theta) 0;...
		0 0 1]*net_disp_orig;

	
end

% Function to integrate up system velocities using a fixed-step method
function [net_disp_orig, cost] = fixed_step_integrator(s,gait,tspan,ConnectionEval,resolution)

	% Duplicate 'resolution' to 'res' if it is a number, or place res at a
	% starting resolution if an automatic convergence method is selected
	% (automatic convergence not yet enabled)
	default_res = 100;
	if isnumeric(resolution)
		res = resolution;
	elseif ischar(resolution) && strcmp(resolution,'autoconverge')
		res = default_res;
	else
		error('Unexpected value for resolution');
	end
	
	% Generate the fixed points from the time span and resolution
	tpoints = linspace(tspan(1),tspan(2),res);
	tsteps = gradient(tpoints);
    
    %Prep interpolation inputs for velocities function
    shape = zeros(size(s.grid.eval));
    shape_gait_def_0 = readGait(gait.phi_def,0);
    actual_size = min(numel(shape),numel(shape_gait_def_0));
    
    samplePoints = {};
    for dim = 1:actual_size
        samplePoints{dim} = [];
    end
    
    for time = tpoints
        shape_gait_def = readGait(gait.phi_def,time);
        shape(1:actual_size) = shape_gait_def(1:actual_size);
        for dim = 1:numel(shape)
            samplePoints{dim}(end+1) = shape(dim); 
        end
    end
    
    indexList = 1:numel(tpoints);
    id = eye(actual_size);
    
    As = cellfun(@(C) -interpn(s.grid.eval{:},C,...
        samplePoints{:},'spline'),s.vecfield.eval.content.Avec,...
        'UniformOutput',false);
    As = celltensorconvert(As);
    
    switch s.costfunction
        case {'pathlength coord','acceleration coord'}
            %In the case of these two cost functions, we only care about
            %the connection field, and the metric is always identity.
            %dM is passed a filler value
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,As{i},id,1),...
                tpoints,indexList,'UniformOutput',false);
        case {'torque','covariant acceleration','power quality'}
            %In the inertial cases, we need to calculate dM, and the metric
            metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
                'UniformOutput',false);
            metrics = celltensorconvert(metrics);
            
            dM_set = {};
            for dim = 1:actual_size
                dM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                    samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.dM{dim},...
                    'UniformOutput',false);
                dM_holder = celltensorconvert(dM_holder);
                dM_set{dim} = dM_holder;
            end
            dMs = {};
            for i = 1:numel(tpoints)
                dMs{i} = {};
                for dim = 1:actual_size
                    dMs{i}{dim} = dM_set{dim}{i};
                end
            end
            
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
                As{i},metrics{i},dMs{i}),tpoints,indexList,'UniformOutput',false);
        otherwise
            %Otherwise, we're not doing inertial so don't need dM, but we
            %do care about the metric and connection
            metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
                'UniformOutput',false);
            metrics = celltensorconvert(metrics);
            
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
                As{i},metrics{i},1),tpoints,indexList,'UniformOutput',false);
    end

	%%%%%%%
	% Integrate cost and displacement into final values
	
	%%
	% Exponential integration for body velocity
	
	% Exponentiate each velocity over the corresponding time step
	expXi = cellfun(@(xi,timestep) se2exp(xi*timestep),xi,num2cell(tsteps),'UniformOutput',false);
	
	% Start off with zero position and displacement
	net_disp_matrix = eye(size(expXi{1}));
	
	% Loop over all the time steps from 1 to n-1, multiplying the
	% transformation into the current displacement
	for i = 1:(length(expXi)-1)
		
		net_disp_matrix = net_disp_matrix * expXi{i};
		
	end
	
	% De-matrixafy the result
	g_theta = atan2(net_disp_matrix(2,1),net_disp_matrix(1,1));
	g_xy = net_disp_matrix(1:2,3);
	
	net_disp_orig = [g_xy;g_theta];
	
	%%
	% Trapezoidal integration for cost
	dcost = [dcost{:}];
	cost = trapz(tpoints,dcost);

end


% Function to evaluate velocity and differential cost at each time for ODE
% solver
function dX = helper_function(t,X,s,gait,ConnectionEval)

	% X is the accrued displacement and cost

	[xi, dcost] = get_velocities(t,s,gait,ConnectionEval);
		
	% Rotate body velocity into world frame
	theta = X(3);
	v = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
		
	% Combine the output
	dX = [v;dcost];
	

end

function expXi = se2exp(xi)

	% Make sure xi is a column
	xi = xi(:);

	% Special case non-rotating motion
	if xi(3) == 0
		
		expXi = [eye(2) xi(1:2); 0 0 1];
		
	else
		
		z_theta = xi(3);
		
		z_xy = 1/z_theta * [sin(z_theta), 1-cos(z_theta); cos(z_theta)-1, sin(z_theta)] * xi(1:2);
		
		expXi = [ [cos(z_theta), -sin(z_theta); sin(z_theta), cos(z_theta)], z_xy;
			0 0 1];
		
	end


end

% function [g_end_orig,g_end_opt, cost_end] = extract_displacement_and_cost(datafile)
% % Extract the displacement and cost data from a sysf_...shchf_....mat file
% 
% % Load the target file
% load(datafile,'p')
% 
% % Prime arrays to hold the net displacement (in original and optimal
% % coordinates) and cost from each shape change in the file. p.G_locus_full is
% % single-level cell array of structures, each of which holds the
% % information for one gait (with all segments concatenated)
% g_end_orig = zeros(numel(p.G_locus_full),3);
% g_end_opt = g_end_orig;
% cost_end = zeros(numel(p.G_locus_full,1)); % If distance metric was not specified, euclidean metric in the parameters was assumed
% 
% % Loop over each shape change
% for i = 1:numel(p.G_locus_full)
% 	
% 	% Extract the final values for the relevant parameters
% 	g_end_orig(i,:) = p.G_locus_full{i}.G(end,:); 
% 	g_end_opt(i,:) = p.G_locus_full{i}.G_opt(end,:); 
% 	cost_end(i) = p.G_locus_full{i}.S(end);
% end
% end

function [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g)
% Calculates the gradient of the shape position, velocity, and acceleration
% with respect to the fourier coefficients.
% Inputs:
%   n: Number of points composing the gait when using shape-space
%       parametrization
%   y: Set of fourier coefficients parametrizing the gait
%   g: Time period of the gait

% Get the fourier frequency, number of shape variables, and number of
% fourier coefficientsinertia_cost_gradient
w = y(end,:);
dim = size(y,2);
num_coeffs = size(y,1);
% Initialize cell array of function handles to hold the partials of the
% shape variables with respect to the fourier coefficients
grad_alpha = cell(num_coeffs,dim);
grad_alphadot = cell(num_coeffs,dim);
grad_alphaddot = cell(num_coeffs,dim);

for i = 1:num_coeffs 
    for j = 1:dim
        % Coefficient a_0 is a lone scalar, so partials with respect to it
        % are zero for the dotted terms and 1 for alpha
        if i == 1
            grad_alpha{i,j} = @(t) [0*(1:j-1), 1, 0*(j+1:dim)];
            grad_alphadot{i,j} = @(t) zeros(1,dim);
            grad_alphaddot{i,j} = @(t) zeros(1,dim);
            continue
        elseif i == num_coeffs % Partial w.r.t. frequency
            grad_alpha{i,j} = @(t) [0*(1:j-1), t*(-y(2,j)*sin(w(j)*t) + y(3,j)*cos(w(j)*t) - ...
                                               2*y(4,j)*sin(2*w(j)*t) + 2*y(5,j)*cos(2*w(j)*t) - ...
                                               3*y(6,j)*sin(3*w(j)*t) + 3*y(7,j)*cos(3*w(j)*t) - ...
                                               4*y(8,j)*sin(4*w(j)*t) + 4*y(9,j)*cos(4*w(j)*t)), ...
                                    0*(j+1:dim)];
            
            grad_alphadot{i,j} = @(t) [0*(1:j-1), -y(2,j)*(w(j)*t*cos(w(j)*t)+sin(w(j)*t)) + y(3,j)*(-w(j)*t*sin(w(j)*t)+cos(w(j)*t)) + ...
                                                  -2*y(4,j)*(2*w(j)*t*cos(2*w(j)*t)+sin(2*w(j)*t)) + 2*y(5,j)*(-2*w(j)*t*sin(2*w(j)*t)+cos(2*w(j)*t)) + ...
                                                  -3*y(6,j)*(3*w(j)*t*cos(3*w(j)*t)+sin(3*w(j)*t)) + 3*y(7,j)*(-3*w(j)*t*sin(3*w(j)*t)+cos(3*w(j)*t)) + ...
                                                  -4*y(8,j)*(4*w(j)*t*cos(4*w(j)*t)+sin(4*w(j)*t)) + 4*y(9,j)*(-4*w(j)*t*sin(4*w(j)*t)+cos(4*w(j)*t)), ...
                                      0*(j+1:dim)];
            
            grad_alphaddot{i,j} = @(t) [0*(1:j-1), -y(2,j)*w(j)*(-t*w(j)*sin(w(j)*t)+2*cos(w(j)*t)) - y(3,j)*w(j)*(t*w(j)*cos(w(j)*t)+2*sin(w(j)*t)) + ...
                                                   -4*y(4,j)*w(j)*(-2*t*w(j)*sin(2*w(j)*t)+2*cos(2*w(j)*t)) - 4*y(5,j)*w(j)*(2*t*w(j)*cos(2*w(j)*t)+2*sin(2*w(j)*t)) + ...
                                                   -9*y(6,j)*w(j)*(-3*t*w(j)*sin(3*w(j)*t)+2*cos(3*w(j)*t)) - 9*y(7,j)*w(j)*(3*t*w(j)*cos(3*w(j)*t)+2*sin(3*w(j)*t)) + ...
                                                   -16*y(8,j)*w(j)*(-4*t*w(j)*sin(4*w(j)*t)+2*cos(4*w(j)*t)) - 16*y(9,j)*w(j)*(4*t*w(j)*cos(4*w(j)*t)+2*sin(4*w(j)*t)), ...
                                       0*(j+1:dim)];
            continue
        end
        % For partial alpha, a_n is associated with cosine and b_n is
        % associated with sine; a_n terms are every second row entry in y
        % with the b_n terms in between
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        % mult comes from the multiplier of the natural frequency for
        % increasing fourier coefficients
        mult = floor(i/2);
        
        grad_alpha{i,j} = @(t) [0*(1:j-1), trig(mult*w(j)*t), 0*(j+1:dim)];
        % For partial alphadot, a_n is associated with sine and b_n is
        % associated with cosine; the a_n terms are every second row entry 
        % in y with the b_n terms in between
        if mod(i,2) == 0
            trig = @sin;
        else
            trig = @cos;
        end
        
        grad_alphadot{i,j} = @(t) [0*(1:j-1), (-1)^(i-1)*mult*w(j)*trig(mult*w(j)*t), 0*(j+1:dim)];
        
        % For partial alphaddot, a_n is associated with cosine and b_n is
        % associated with sine
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        
        grad_alphaddot{i,j} = @(t) [0*(1:j-1), -mult^2*w(j)^2*trig(mult*w(j)*t), 0*(j+1:dim)];
    end
end
end

function cost_grad = inertia_cost_gradient(s,n,y,g,gait,EvaluationMethod)
% Calculates the gradient of cost for inertial systems.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   EvaluationMethod: String representing how the gradient of cost should
%   be calculated; provide 'discrete' to evaluate at 100 discrete values or
%   'ode45' if you would like the gradient of cost to be integrated using
%   ode45. Other values will result in error.

    % Contribution to gradient from the movement of each point due to
    % change in fourier coefficients
    [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);

    cost_grad = zeros(size(grad_alpha));
    tspan = [0 g];
    if strcmpi(EvaluationMethod,'discrete')
        num_pts = 100;
        t_pts = linspace(0,g,num_pts);
        del_t = t_pts(2) - t_pts(1);
        
        %Prep interpolation inputs for gradient calcs
        shape = zeros(size(s.grid.eval));
        shape_gait_def_0 = readGait(gait.phi_def,0);
        actual_size = min(numel(shape),numel(shape_gait_def_0));
        
        samplePoints = {};
        for dim = 1:actual_size
            samplePoints{dim} = [];
        end
        
        for time = t_pts
            shape_gait_def = readGait(gait.phi_def,time);
            shape(1:actual_size) = shape_gait_def(1:actual_size);
            for dim = 1:numel(shape)
                samplePoints{dim}(end+1) = shape(dim);
            end
        end
        
        %Batch interpolate the metric at each point along the gait
        metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
            samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
            'UniformOutput',false);
        metrics = celltensorconvert(metrics);
        
        dM_set = {};
        for dim = 1:actual_size
            dM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.dM{dim},...
                'UniformOutput',false);
            dM_holder = celltensorconvert(dM_holder);
            dM_set{dim} = dM_holder;
        end
        dMs = {};
        for i = 1:numel(t_pts)
            dMs{i} = {};
            for dim = 1:actual_size
                dMs{i}{dim} = dM_set{dim}{i};
            end
        end
        
        empty_ddM = cell(size(s.coriolisfield.coriolis_eval.content.ddM));
        ddM_set = empty_ddM;
        for dim = 1:numel(ddM_set)
            ddM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.ddM{dim},...
                'UniformOutput',false);
            ddM_holder = celltensorconvert(ddM_holder);
            ddM_set{dim} = ddM_holder;
        end
        ddMs = {};
        for i = 1:numel(t_pts)
            ddMs{i} = empty_ddM;
            for dim = 1:numel(empty_ddM)
                ddMs{i}{dim} = ddM_set{dim}{i};
            end
        end
                
        switch s.costfunction
            case 'torque'
                for k = 1:length(t_pts)
                    del_cost = inertia_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'covariant acceleration'
                for k = 1:length(t_pts)
                    del_cost = acceleration_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'acceleration coord'
                for k = 1:length(t_pts)
                    del_cost = accelerationcoord_gradient_helper(t_pts(k),[],s,gait,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'power quality'
                for k = 1:length(t_pts)
                    del_cost = powerquality_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
        end
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    elseif strcmpi(EvaluationMethod,'ode45')
        switch s.costfunction
            case 'torque'
                sol = ode45(@(t,y) inertia_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'covariant acceleration'
                sol = ode45(@(t,y) acceleration_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'acceleration coord'
                sol = ode45(@(t,y) accelerationcoord_gradient_helper(t,y,s,gait,grad_alphaddot),tspan,cost_grad);
            case 'power quality'
                sol = ode45(@(t,y) powerquality_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
        end
        % Extract the final motion
        cost_grad = reshape(deval(sol,tspan(end)),size(cost_grad));
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    else
        error('Untenable option provided for EvaluationMethod!')
    end
end

function validate_shape_gradient(n,y,g,grad_alphaddot,grad_alphadot,grad_alpha) %#ok<DEFNU>
% Function that helps validate that the gradient of shape position,
% velocity, and acceleration are correctly calculated. Difference between
% the input gradients and calculation-verified gradients are printed to the
% terminal. Should be very close to zero.
% Inputs:
%   n: Number of points at which gait should be evaluated in shape space
%   y: Fourier coefficients that parametrize the gait
%   g: Time period over which gait is executed
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

    % Select 10 random times at which to evaluate the gradient of shape
    t = g*rand(1,10);
    w1 = y(end,1); % Frequency of Fourier transform
    w2 = y(end,2);
    
    % Evaluate the gradient at each random time against a hard-coded
    % gradient calculation
    for i = 1:length(t)
        grad_alpha_eval = cellfun(@(C) C(t(i)), grad_alpha, 'UniformOutput', false);
        grad_alpha1_eval = cell2mat(grad_alpha_eval(:,1));
        grad_alpha2_eval = cell2mat(grad_alpha_eval(:,2));
        grad_alpha1_calc = [1,0;cos(w1*t(i)),0;sin(w1*t(i)),0;cos(2*w1*t(i)),0;sin(2*w1*t(i)),0;cos(3*w1*t(i)),0;sin(3*w1*t(i)),0;cos(4*w1*t(i)),0;sin(4*w1*t(i)),0;0,0];
        grad_alpha2_calc = [0,1;0,cos(w2*t(i));0,sin(w2*t(i));0,cos(2*w2*t(i));0,sin(2*w2*t(i));0,cos(3*w2*t(i));0,sin(3*w2*t(i));0,cos(4*w2*t(i));0,sin(4*w2*t(i));0,0];
        grad_alpha1_err = grad_alpha1_eval - grad_alpha1_calc
        grad_alpha2_err = grad_alpha2_eval - grad_alpha2_calc
        
        grad_alphadot_eval = cellfun(@(C) C(t(i)), grad_alphadot, 'UniformOutput', false);
        grad_alphadot1_eval = cell2mat(grad_alphadot_eval(:,1));
        grad_alphadot2_eval = cell2mat(grad_alphadot_eval(:,2));
        grad_alphadot1_calc = [0,0;-w1*sin(w1*t(i)),0;w1*cos(w1*t(i)),0;-2*w1*sin(2*w1*t(i)),0;2*w1*cos(2*w1*t(i)),0;-3*w1*sin(3*w1*t(i)),0;3*w1*cos(3*w1*t(i)),0;-4*w1*sin(4*w1*t(i)),0;4*w1*cos(4*w1*t(i)),0;0,0];
        grad_alphadot2_calc = [0,0;0,-w2*sin(w2*t(i));0,w2*cos(w2*t(i));0,-2*w2*sin(2*w2*t(i));0,2*w2*cos(2*w2*t(i));0,-3*w2*sin(3*w2*t(i));0,3*w2*cos(3*w2*t(i));0,-4*w2*sin(4*w2*t(i));0,4*w2*cos(4*w2*t(i));0,0];
        grad_alphadot1_err = grad_alphadot1_eval - grad_alphadot1_calc;
        grad_alphadot2_err = grad_alphadot2_eval - grad_alphadot2_calc;
        
        grad_alphaddot_eval = cellfun(@(C) C(t(i)), grad_alphaddot, 'UniformOutput', false);
        grad_alphaddot1_eval = cell2mat(grad_alphaddot_eval(:,1));
        grad_alphaddot2_eval = cell2mat(grad_alphaddot_eval(:,2));
        grad_alphaddot1_calc = [0,0;-w1^2*cos(w1*t(i)),0;-w1^2*sin(w1*t(i)),0;-4*w1^2*cos(2*w1*t(i)),0;-4*w1^2*sin(2*w1*t(i)),0;-9*w1^2*cos(3*w1*t(i)),0;-9*w1^2*sin(3*w1*t(i)),0;-16*w1^2*cos(4*w1*t(i)),0;-16*w1^2*sin(4*w1*t(i)),0;0,0];
        grad_alphaddot2_calc = [0,0;0,-w2^2*cos(w2*t(i));0,-w2^2*sin(w2*t(i));0,-4*w2^2*cos(2*w2*t(i));0,-4*w2^2*sin(2*w2*t(i));0,-9*w2^2*cos(3*w2*t(i));0,-9*w2^2*sin(3*w2*t(i));0,-16*w2^2*cos(4*w2*t(i));0,-16*w2^2*sin(4*w2*t(i));0,0];
        grad_alphaddot1_err = grad_alphaddot1_eval - grad_alphaddot1_calc;
        grad_alphaddot2_err = grad_alphaddot2_eval - grad_alphaddot2_calc;
    end
end

function validate_cost_gradient(s,n,y,g,p)
% Function to help verify that the cost gradient calculated by
% inertia_gradient_helper is valid. Values returned by function and
% individually calculated costs are printed to the terminal along with the
% difference between the two methods. Should be relatively small for full
% range of fourier coefficients.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   p: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t

% Set the delta in the fourier coefficients between individual cost
% evaluations
fourier_delta = 0.001;
[grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);
% Go through each of the fourier coefficients and test how applying
% fourier_delta to each of them matches with the cost gradient calculated
% by inertia_gradient_helper
for fourier_test = 1:numel(y)
    % Perturb the fourier coefficient to be tested by fourier_delta and
    % create a parametrization for calculating the point in the gait at a
    % time t
    y2 = y;
    y2(fourier_test) = y2(fourier_test) + fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p2 = makeGait(y2);
            
    % Perturb again in the opposite direction
    y2 = y;
    y2(fourier_test) = y2(fourier_test) - fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p3 = makeGait(y2);
     
    % Get the shape and shape derivative at a random time for each gait
    t = g*rand(1);
	shape = readGait(p3.phi_def,t);
	shapelist = num2cell(shape);
	dshape = readGait(p3.dphi_def,t);
    ddshape = readGait(p3.ddphi_def,t);
    shape_delta = readGait(p2.phi_def,t);
	shapelist_delta = num2cell(shape_delta);
	dshape_delta = readGait(p2.dphi_def,t);
    ddshape_delta = readGait(p2.ddphi_def,t);
    
    % Get mass matrices at both locations
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    M_delta = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist_delta{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    
    % Get partial mass matrices at both locations
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    dM_alphadalpha_delta = calc_partial_mass(s,shapelist_delta);
    
    % Calculate the cost using the two different gaits
    cost = torque_cost(M,dM_alphadalpha,shape,dshape,ddshape);
    cost_delta = torque_cost(M_delta,dM_alphadalpha_delta,shape_delta,dshape_delta,ddshape_delta);
    % Calculate what the gradient of cost is for this particular point
    cost_grad = inertia_gradient_helper(t,[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
    cost_grad = reshape(cost_grad,size(y));
    cost_grad_rel = cost_grad(fourier_test)
    cost_grad_calc = (cost_delta-cost)/(2*fourier_delta)
    
    % Find what the difference is between cost_grad and the costs evaluated
    % at distance fourier_delta
    err = cost_grad_rel - cost_grad_calc
end

end

function del_cost = inertia_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of cost for inertia systems.
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
    shape = zeros(size(s.grid.eval));
    dshape = zeros(size(s.grid.eval));
    ddshape = zeros(size(s.grid.eval));
    
	shape_gait_def = readGait(gait.phi_def,t);
	dshape_gait_def = readGait(gait.dphi_def,t);
    ddshape_gait_def = readGait(gait.ddphi_def,t);
    
    actual_size = min(numel(shape),numel(shape_gait_def));
    shape(1:actual_size) = shape_gait_def(1:actual_size);
    dshape(1:actual_size) = dshape_gait_def(1:actual_size);
    ddshape(1:actual_size) = ddshape_gait_def(1:actual_size);

    shapelist = num2cell(shape);
    
    metricgrad = getMetricGrad(s,shape,dM,grad_alpha_eval);
    
    M = metric;
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM,shape,dshape);
    tau = M*ddshape(:) + C;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM{j}*del_shape(j);
        end
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM{j}*del_dshape(:);
        end
        
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        del_cost(i) = del_tau(:)'*tau(:)...
                    + tau(:)'*del_tau(:);
    end
    del_cost = del_cost(:);
end

function del_cost = acceleration_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of covariant acceleration cost
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
    shape = zeros(size(s.grid.eval));
    dshape = zeros(size(s.grid.eval));
    ddshape = zeros(size(s.grid.eval));
    
	shape_gait_def = readGait(gait.phi_def,t);
	dshape_gait_def = readGait(gait.dphi_def,t);
    ddshape_gait_def = readGait(gait.ddphi_def,t);
    
    actual_size = min(numel(shape),numel(shape_gait_def));
    shape(1:actual_size) = shape_gait_def(1:actual_size);
    dshape(1:actual_size) = dshape_gait_def(1:actual_size);
    ddshape(1:actual_size) = ddshape_gait_def(1:actual_size);
            
	shapelist = num2cell(shape);
    
    metricgrad = getMetricGrad(s,shape,dM,grad_alpha_eval);
    
    M = metric;
    
    % Using inverse enough that it's worth it to calc once upfront
    M_inv = inv(M);
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM,shape,dshape);
    tau = M*ddshape(:) + C;
    
    %Covariant acceleration calculation
    cov_acc = M_inv*tau;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM{j}*del_shape(j);
        end
        
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        
        %Gradient of inverse of mass for covariant acc calculation
        Minv_grad = zeros(length(shapelist));
        for j = 1:length(shapelist)
            %Formula from matrix cookbook
            %d(M^-1)=-M^-1*dM*M^-1
            dM_alphainv_dalpha = -M_inv*dM{j}*M_inv;
            Minv_grad = Minv_grad + dM_alphainv_dalpha*del_shape(j);
        end
        
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM{j}*del_dshape(:);
        end
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        % Gradient of covariant acc
        del_cov_acc = Minv_grad*tau(:)+M_inv*del_tau(:);
        del_cost(i) = del_cov_acc'*metric*cov_acc...
                    + cov_acc'*metricgrad{i}*cov_acc...
                    + cov_acc'*metric*del_cov_acc;
    end
    del_cost = del_cost(:);
end

function del_cost = accelerationcoord_gradient_helper(t,X,s,gait,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of shape-space acceleration
% cost
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alphaddot_eval));
    % Get the shape derivative at the current time
    
    ddshape_def = readGait(gait.ddphi_def,t);
    actual_size = min(numel(s.grid.eval),numel(ddshape_def));
    ddshape = zeros(size(s.grid.eval));
    
    ddshape(1:actual_size) = ddshape_def(1:actual_size);
   
    % Regular cost calculation
    cost = sqrt(ddshape(:)'*ddshape(:));
    
    for i = 1:numel(grad_alphaddot_eval)
        % Partial of shape acceleration with respect to fourier coefficient i
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of shape-space acceleration cost
        del_cost(i) = del_ddshape(:)'*ddshape(:)+ddshape(:)'*del_ddshape(:);
    end
    del_cost = del_cost(:);
end

function del_cost = powerquality_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of cost for inertia systems.
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
	shape = readGait(gait.phi_def,t);
	shapelist = num2cell(shape);
	dshape = readGait(gait.dphi_def,t);
    ddshape = readGait(gait.ddphi_def,t);
    
    % Get mass and partial mass matrices
    M = metric;
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM,shape,dshape);
    tau = M*ddshape(:) + C;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM{j}*del_shape(j);
        end
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM{j}*del_dshape(:);
        end
        
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        % Gradient of (P1+P2)^2, aka gradient of v'*tau*v'*tau where v
        % is the shape velocity
        qual_1 = del_dshape(:)'*tau*dshape(:)'*tau + ...
                 dshape(:)'*del_tau*dshape(:)'*tau + ...
                 dshape(:)'*tau*del_dshape(:)'*tau + ...
                 dshape(:)'*tau*dshape(:)'*del_tau;
        % Gradient of P1^2 + P2^2
        qual_2 = (2*del_dshape(:).*dshape(:))'*(tau.*tau) + ...
                 (dshape(:).*dshape(:))'*(2*del_tau.*tau);
        % Gradient of power quality cost wrt this fourier coeff
        del_cost(i) = qual_1 - qual_2;
    end
    del_cost = del_cost(:);
end

function metricgrad = getMetricGrad(s,shape,dM,grad_alpha)
%Returns metric at a shape position, and gradient of metric w.r.t. fourier
%coefficients

%s - structure containing metric function
%shape - array of shape values
%dM - derivative of matrix with respect to shape variables
%grad_alpha - gradient of shape values w.r.t. fourier coefficients

    n_dim = numel(s.grid.eval);
    actual_size = min(size(shape,2),n_dim);
    
    metricgrad = repmat({zeros(actual_size)},size(grad_alpha));
    for j = 1:actual_size
        
        for k = 1:numel(grad_alpha)
            metricgrad{k} = metricgrad{k}+ dM{j}*grad_alpha{k}(j);
        end
    end

end