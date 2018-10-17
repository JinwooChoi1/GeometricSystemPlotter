function output = sysf_four_link_lowRe(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Viscous Swimmer: 4-link'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/N_link_chain.m';...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_discrete.m';...
				'Utilities/LowRE_toolbox/LowRE_local_connection_discrete.m'});

		case 'initialize'

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            
            
            %%%
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. The code below uses properties of cell
            % arrays to automatically match the dimensionality of the grid
            % with the number of shape basis functions in use
            s.visual.grid = cell(numel(s.geometry.linklengths)-1,1);
            [s.visual.grid{:}] = ndgrid([-1  0  1]);

            
            %%%
            %%%%%%
            % Define system physics
            s.physics.drag_ratio = 2;
            s.physics.drag_coefficient = 1;
           
 
            %Functional Local connection and dissipation metric

            s.A = @(alpha1,alpha2,alpha3) LowRE_local_connection( ...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,alpha3]);                        % Joint angles
            
            s.metric = @(alpha1,alpha2,alpha3) LowRE_dissipation_metric(...
                        s.geometry,...                           % Geometry of body
                        s.physics,...                            % Physics properties
                        [alpha1,alpha2,alpha3]);                        % Joint angles

                    
			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1,-1,1]*2;

			%densities for various operations
			s.density.vector = [10 10 10]; %density to display vector field
			s.density.scalar = [21 21 21]; %density to display scalar functions
			s.density.eval = [21 21 21];   %density for function evaluations
            s.density.metric_eval = [11 11 11]; %density for metric evaluation
            s.density.finite_element=11;

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*1;
			s.tic_locs.y = [-1 0 1]*1;


			%%%%
			%Save the system properties
			output = s;


	end

end
