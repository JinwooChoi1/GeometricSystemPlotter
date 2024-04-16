function vel_plot_from_fourier()
    %Store handle to current axis so we can do edits
    newfigurehandle = gca;

    %Prep figure from sysplotter
    axis tight;
    hold on;
    %colormap(cmap)

    %Set arrow prongs to a good thickness
    linenum = size(newfigurehandle.Children,1) - 1;  
    
    lineidxset = [];
    arrowidxset = [];
    for i = 1:linenum
        if length(newfigurehandle.Children(i).XData) == 50
            lineidxset = [lineidxset i];
        elseif length(newfigurehandle.Children(i).XData) == 2
            arrowidxset = [arrowidxset i];
        end
    end

    gaitnum = length(lineidxset);
    
    XData = cell(gaitnum,1);
    YData = cell(gaitnum,1);
    for gaitidx = 1:gaitnum

        %Get shape data over gait
        XData{gaitidx} = newfigurehandle.Children(lineidxset(gaitidx)).XData;
        YData{gaitidx} = newfigurehandle.Children(lineidxset(gaitidx)).YData;

    end

    %Remove the line of constant width
    delete(newfigurehandle.Children([lineidxset arrowidxset]));

    for gaitidx = 1:gaitnum
        %Calculate distances between points on gait
        dXData = diff(XData{gaitidx});
        dYData = diff(YData{gaitidx});

        %Calculate length of each line segment as fraction of whole
        %gait
        percentages = sqrt(dXData.^2+dYData.^2)/...
            sum(sqrt(dXData.^2+dYData.^2));

        %Calculate min and max of these ratios
        minp = min(percentages);
        maxp = max(percentages);

        %Set linewidths so that slower motion is a thicker line
        linewidths = interp1([minp,maxp],[7,3],percentages);

        %Redraw gait with varied linewidths
        for i = 1:49
            plot(XData{gaitidx}(i:i+1),YData{gaitidx}(i:i+1),'-r','LineWidth',linewidths(i));
        end

        %Draw circles under the joints between line segments to make it
        %less obvious what we did
        scatter(XData{gaitidx}(1:49),YData{gaitidx}(1:49),pi*(linewidths/2).^2,'r','filled');
    end
end