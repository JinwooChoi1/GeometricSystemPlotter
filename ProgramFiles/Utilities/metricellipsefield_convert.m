function h = metricellipsefield_convert(x,y,M,style,convert,stretch,varargin)
    %Plot an ellipse field based on the singular values of a metric
    %tensor M. M should be specified as a cell array with with the same
    %dimensions as x and y, and each cell containing a 2x2 matrix containing
    %its value at the corresponding x,y location


    %%%%%%%%%%%%%
    % Construct the ellipses for each location

    % Make the circle primitive
    th = linspace(0,2*pi,50);
    xc = cos(th);
    yc = sin(th);

    % Replicate the circle primitive into a cell array at each x,y location
    circles = repmat({[xc;yc]},size(x));

    % Pull out the stretch name
    stretchnames = {'stretch','surface'};
    if stretch
        stretchname = stretchnames{stretch};
    end


    if stretch
        J = arrayfun(@(u,v) convert.(stretchname).jacobian(u,v),x,y,'UniformOutput',false);
    end

    if strcmpi(varargin{3},'M')
        [u,s,~] = cellfun(@(m) svd(inv(m)),M,'UniformOutput',false);
    elseif strcmpi(varargin{3},'Co')
        if stretch
            % Calculate the svd of the cometric tensor for metric stertch
            [u,s,~] = cellfun(@(m,j) svd((j.')^-1*eye(size(m))*j^-1),M,J,'UniformOutput',false);
        else
            % Calculate the svd of the cometric tensor
            [u,s,~] = cellfun(@(m) svd(eye(size(m))),M,'UniformOutput',false);
        end
    elseif strcmpi(varargin{3},'Dual')
        if stretch
            % Calculate the svd of the cometric tensor for metric stertch
            J = cellfun(@(u)norm(u)*u^-1,J,'UniformOutput',false);
        end
        % Calculate the svd of the cometric tensor
        [u,s,~] = cellfun(@(m) svd(m),M,'UniformOutput',false);
    else
        error("metric field is not chosen")
    end

    % Apply the transform corresponding to M to the circles
    circles = cellfun(@(u,s,c)(u*sqrt(s))*c,u,s,circles,'UniformOutput',false);

    switch style

    	% The tissot indicatrix, showing linear stretch
    	case 'tissot'

    		plot_options = varargin{1};

        	% the tissot indicatrix, with major and minor axis crosses added
    	case 'tissot-cross'

    		% Create the crosses
    		crosses = cellfun(@(u,s)...
    			u*[-sqrt(s(1,1)) sqrt(s(1,1)) NaN 0 0;...[u(1,1)*s(1,1) -u(1,1)*s(1,1);...% NaN u(1,2)*s(2,2) u(1,2)*s(2,2);...
    			0 0 NaN -sqrt(s(2,2)) sqrt(s(2,2))...u(2,1)*s(1,1) -u(2,1)*s(1,1);... NaN u(2,2)*s(2,2) -u(2,2)*s(2,2)
    			]...
                ,u,s,'UniformOutput',false);

    		plot_options = varargin{1};
    		plot_options_crosses = varargin{2};
    end

    % Recenter the circles and stretch them if there is a conversion transform
    if stretch
        % For metric ellipse field, take the Jacobian.
        % For dualmetric ellipse field, take the inverse jacobian but J is
        % already inversed.
        if strcmpi(varargin{3},'M') || strcmpi(varargin{3},'Dual')  
            % Boost Jacobian dimension if needed
            if stretch == 1
                J = cellfun(@(u)[u;0 0],J,'UniformOutput',false);
            end

            circles = cellfun(@(u,v)u*v,J,circles,'UniformOutput',false);
            if exist('crosses','var')
                crosses = cellfun(@(u,v)u*v,J,crosses,'UniformOutput',false);
            end

        % For cometric ellipse field, take the identity matrix because the
        % ellipse already match with the isometric embedding surface.
        elseif strcmpi(varargin{3},'Co')
            circles = cellfun(@(u)[eye(2,2); 0 0]*u,circles,'UniformOutput',false);
            if exist('crosses','var')
                crosses = cellfun(@(u) [eye(2,2); 0 0]*u,crosses,'UniformOutput',false);
            end
        end
    end
    
    if stretch
        [x, y, z] = convert.(stretchname).old_to_new_points(x,y);
    end

    % Find the greatest x or y range in the circles
    xrange = cellfun(@(u)range(u(1,:)),circles);
    yrange = cellfun(@(u)range(u(2,:)),circles);

    xrange = max(xrange(:));
    yrange = max(yrange(:));

    % Find the smallest spacing in each direction
    xspacing = min(diff(x(:,1)));
    yspacing = min(diff(y(1,:)));

    % Set the percentage of the spacing that the largest element should occupy
    max_fill = 0.6;

    % Determine if x or y fitting is the limiting factor for the plot and set
    % the scaling accordingly
    scale_factor = min(xspacing/xrange, yspacing/yrange)*max_fill;

    % Multiply all the circles by the scale factor
    circles = cellfun(@(u)u*scale_factor,circles,'UniformOutput',false);

    if exist('crosses','var')
    	crosses = cellfun(@(u)u*scale_factor,crosses,'UniformOutput',false);
    end
    
    circles = cellfun(@(u,v,w) [u(1,:)+v;u(2,:)+w;zeros(size(u(1,:)))],circles,num2cell(x),num2cell(y),'UniformOutput',false);
    if exist('crosses','var')
        crosses = cellfun(@(u,v,w) [u(1,:)+v;u(2,:)+w;zeros(size(u(1,:)))],crosses,num2cell(x),num2cell(y),'UniformOutput',false);
    end





    %%%%%%%%%%%%%%%%%
    % Fill in default values for plot options

    % Ensure plot options exists
    if ~exist('plot_options','var')

    	plot_options = {};

    end


    % Color
    if ~any(strcmpi('edgecolor',plot_options(1:2:end)))

    	plot_options = [plot_options,{'EdgeColor','k'}];

    end
    if ~any(strcmpi('facecolor',plot_options(1:2:end)))

    	plot_options = [plot_options,{'FaceColor','none'}];

    end

    % Parent
    if isempty(plot_options) || ~any(strmatch('parent',lower(plot_options(1:2:end))))

    	f = figure;
    	ax = axes('Parent',f);
    	plot_options = [plot_options,{'Parent',ax}];

    end

    %%%%%%%%%%%%%%%%%%
    % Make the ellipses

    switch style

    	case 'tissot'

    		h = cellfun(@(u)patch('XData',u(1,:),'YData',u(2,:),'ZData',u(3,:),plot_options{:}),circles,'UniformOutput',false);

    	case 'tissot-cross'

    		%h_cross = cellfun(@(u)patch('XData',u(1,:),'YData',u(2,:),plot_options_crosses{:}),crosses,'UniformOutput',false);


            h_cross=cell(size(crosses));
            for i =1:numel(h_cross)
                h_cross{i} = line('xdata',crosses{i}(1,:),'ydata',crosses{i}(2,:),'zdata',crosses{i}(3,:),plot_options_crosses{:});
            end

    		h_ellipse = cellfun(@(u)patch('XData',u(1,:),'YData',u(2,:),'ZData',u(3,:),plot_options{:}),circles,'UniformOutput',false);

    end


        %     switch stretch
        %
        %         case 1 %metric stretch
        %             % Multiply all the circles by their respective Jacobians
        %             J = arrayfun(@(u,v) convert.(stretchname).jacobian(u,v),x,y,'UniformOutput',false);
        %             circles = cellfun(@(u,v)u*v,J,circles,'UniformOutput',false);
        %             if exist('crosses','var')
        %                 crosses = cellfun(@(u,v)u*v,J,crosses,'UniformOutput',false);
        %             end
        %
        %             [x_f, y_f] = convert.(stretchname).old_to_new_points(x,y);
        %
        %             circles = cellfun(@(u,v,w) [u(1,:)+v;u(2,:)+w],circles,num2cell(x_f),num2cell(y_f),'UniformOutput',false);
        %             if exist('crosses','var')
        %                 crosses = cellfun(@(u,v,w) [u(1,:)+v;u(2,:)+w],crosses,num2cell(x_f),num2cell(y_f),'UniformOutput',false);
        %             end
        %
        %         case 2 %metric surface