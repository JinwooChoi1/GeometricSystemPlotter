function set_tics_shapespace(ax,s)
% Place the tic marks

    % Set the tic fontsize. 
    % For Mac, 30pt seems to be suitable.
    % For Window, 20pt.
    if ismac        
        font_size = 30;
    else
        font_size = 20;
    end
    set(ax,'FontSize',font_size,'FontName','Helvetica')
		
    set(ax,'XTick',s.tic_locs.x,'YTick',s.tic_locs.y)

end