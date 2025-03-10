function [f,l] = plot_from_fourier(y,nfparam,dimension,figurenum,options)
    arguments
        y double
        nfparam (1,1) double
        dimension (1,1) double
        figurenum (1,1) double
        options.LineWidth (1,1) double = 3
        options.Color
        options.LineStyle
        options.Stretch
    end

    f = figure(figurenum);

    if (size(y,1) == 1)
        y= y.';
    end

    if size(y,1) == dimension*nfparam
        y = reshape(y,[nfparam dimension]);
    elseif size(y,1) == dimension*(nfparam - 1)
        y = reshape(y,[nfparam-1 dimension]);
        y = [y; 2*pi*ones(1,dimension)];
    end

    y1=path_from_fourier(y,99,dimension);
    hold on;
    axis equal;

    if isfield(options,'Stretch')

        if size(y1,2) == 2
            [new_y1(:,1),new_y1(:,2)] = options.Stretch.old_to_new_points(y1(:,1),y1(:,2));
            y1=new_y1;
        else
            error("This function does not supprot 3d metric stretch.")
        end
    end

    if size(y1,2) == 2
        l = line(y1(:,1),y1(:,2));
    elseif size(y1,2) == 3
        l = line(y1(:,1),y1(:,2),y1(:,3));
        view(3);
    end

    if isfield(options,'Color')
        set(l,'Color',options.Color);
    end

    set(l,'LineWidth',options.LineWidth);

    if isfield(options,'LineStyle')
        set(l,'LineStyle',options.LineStyle);
    end

    if exist("l","var")
        l.XData = y1(:,1);
        l.YData = y1(:,2);
    end
end

