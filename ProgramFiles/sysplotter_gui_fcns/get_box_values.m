%Get box names and values from the working column    
function [checkbox_names,...
            checkbox_active,...
            checkbox_values,...
            checkbox_enabled,...
            plot_types,...
            plot_subtypes,...
            plot_style,...
            plot_coordinates]...
	= get_box_values(source_number_text,handles)

    %%%%
    %build the list of plot boxes for the column
    
    %classes of plot
    plot_types = {'vfield',...
                  'CCF',...
                  'mfield',...
                  'ShapeSpace',...
                  'disp',...
                  'beta',...
                  'dbeta',...
                  'xy'};
    
    %plot or subplot for each coordinate
    plot_style = {'plot',...                % Vfield
                  'plot',...                % CCF
                  'single',...              % mfield
                  'single',...              % ShapeSpace
                  'subplot',...             % disp
                  'plot',...                % beta
                  'plot',...                % dbeta
                  'single'};                % xy
    
    %text for the different components
    plot_subtypes = { {'X','Y','T'},...               % Vfield
                      {'X','Y','T'},...               % CCF
                      {'M','Co','Dual'},...           % mfield
                      {''},...                        % shapespace
                      {'X','Y','T'},...               % disp
                      {'X','Y','T'},...               % beta
                      {'X','Y','T'},...               % dbeta
                      {'traj','net','BVI','cBVI'}};   % xy
	
	%Raw or optimized versions
	plot_coordinates = { {'','opt'},...     % Vfield
                         {'','opt'},...     % CCF
                         {''},...           % mfield
                         {'','opt'},...           % ShapeSpace
                         {'','opt'},...     % disp
                         {''},...     % beta
                         {''},...     % dbeta
                         {'','opt'}};     % xy;
	
    
    % Iterate through the plot types, getting names, values, and enabled
    % state
    checkbox_names = cell(size(plot_types));
    checkbox_values = checkbox_names;
    checkbox_enabled = checkbox_names;
    checkbox_active = checkbox_names;
    
    for idx = 1:numel(plot_types)
        
        % Iterate through the plot coordinates
        checkbox_names{idx} = cell(size(plot_coordinates{idx}));
        for idx2 = 1:numel(plot_coordinates{idx})
            
            % Iterate through plot subtypes
            checkbox_names{idx}{idx2} = cell(size(plot_subtypes{idx}));
            for idx3 = 1:numel(plot_subtypes{idx})
                
                % Build a name for the checkbox
                checkbox_names{idx}{idx2}{idx3} = ...
                    [plot_types{idx}...
                     plot_subtypes{idx}{idx3}...
                     plot_coordinates{idx}{idx2}...
                     'checkbox'...
                      source_number_text];
                  
                % Get the checkbox's value
                checkbox_values{idx}{idx2}{idx3} = ...
                    get(handles.(checkbox_names{idx}{idx2}{idx3}),'Value');

                % Get the checkbox's enabled state
                checkbox_enabled{idx}{idx2}{idx3} = ...
                    get(handles.(checkbox_names{idx}{idx2}{idx3}),'Enable');
                
                % Active boxes are both checked and enabled
                checkbox_active{idx}{idx2}{idx3} = ...
                   checkbox_values{idx}{idx2}{idx3} && onoff(checkbox_enabled{idx}{idx2}{idx3});
 
            end
            
        end
        
    end
    
    % Count how many checkboxes are enabled
    checkbox_active_flat = [checkbox_active{:}];
    checkbox_active_flat = [checkbox_active_flat{:}];
    num_checks = sum([checkbox_active_flat{:}]);
    
    if num_checks == 0
        warning('No checkboxes selected.')
        %             return
    end
        
    
% 	%merge plot subtypes
% 	for k = 1:length(plot_subtypes)
% 		for i = 1:length(connection_types)
% 			for j = 1:length(plot_subtypes{k})
% 				merged_plot_subtypes{k}{(i-1)*length(plot_subtypes{k})+j}...
% 					= [plot_subtypes{k}{j} connection_types{i}];
% 			end
% 		end
% 	end
% 	
% 	% Trim out redundant names from the beta and dbeta boxes, which don't
% 	% have opt versions
% 	merged_plot_subtypes{strcmp('beta',plot_types)}...
% 		(cellfun(@(x) ~isempty(x),(strfind(...
% 		merged_plot_subtypes{strcmp('beta',plot_types)},'opt')))) = [];
% 	
% 	merged_plot_subtypes{strcmp('dbeta',plot_types)}...
% 		(cellfun(@(x) ~isempty(x),(strfind(...
% 		merged_plot_subtypes{strcmp('dbeta',plot_types)},'opt')))) = [];
%     
% 	% Trim out redundant names from the xy plot boxes (which handle
% 	% optimized and not optimized more separately
% 	merged_plot_subtypes{strcmp('xy',plot_types)}...
% 		(cellfun(@(x) ~isempty(x),(strfind(...
% 		merged_plot_subtypes{strcmp('xy',plot_types)},'opt')))) = [];
% 	
% 	merged_plot_subtypes{strcmp('xyopt',plot_types)}...
% 		(cellfun(@(x) ~isempty(x),(strfind(...
% 		merged_plot_subtypes{strcmp('xyopt',plot_types)},'opt')))) = [];
% 
%     merged_plot_subtypes{strcmp('mfield',plot_types)}...
% 		(cellfun(@(x) ~isempty(x),(strfind(...
% 		merged_plot_subtypes{strcmp('mfield',plot_types)},'opt')))) = [];
% 
% 	
% 	%prime the arrays of box names and values
%     box_names = cell(1,length(plot_types));%(2*length(plot_subtypes)+1)
%     box_values = cell(size(box_names));
%     box_enabled = cell(size(box_names));
%     
%     %build the list
%     for i = 1:length(plot_types)
% 		
% 		box_names{i} = cell(1+length(merged_plot_subtypes{i}),1);
% 		box_values{i} = zeros(size(box_names{i}));
% 		box_enabled{i} = zeros(size(box_names{i}));                
% 
% 		%category name
%         box_names{i}{1} = plot_types{i};
%         
%         %category value
%         box_values{i}(1) = ...
%                 get(handles.([plot_types{i} 'checkbox' source_number_text]),'Value');
%             
%         %category enabled state
%         box_enabled{i}(1) = ...
% 				strcmp('on',get(handles.([plot_types{i} 'checkbox' source_number_text]),'Enable'));
%             
% 
%         %specific plot names
%         for j = 1:length(merged_plot_subtypes{i})
% 			
% 			%name
% 			box_names{i}{1+j} = ...
% 				[plot_types{i} plot_subtypes{j}];
% 
% 			%value
% 			box_values{i}(1+j) = ...
% 				get(handles.([plot_types{i} merged_plot_subtypes{i}{j} 'checkbox' source_number_text]),'Value');
% 
% 			%enabled state
% 			box_enabled{i}(1+j) = ...
% 				strcmp('on',get(handles.([plot_types{i} merged_plot_subtypes{i}{j} 'checkbox' source_number_text]),'Enable'));
% 			
%         end
%         
%     end
    
%     %Active boxes are both checked and enabled
%     box_active = cellfun(@(x,y) x & y,box_values , box_enabled,'UniformOutput',false);
% 	
%     num_checks = 0; 
%     for i = 1:numel(box_values)
%         
%         num_checks = num_checks + sum(box_values{1,i}(2:end));
% 
%     end
% 
%     if num_checks == 0
%         error('Warning: No checkboxes selected.')
%         %             return
%     end
    
    
end