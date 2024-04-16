function f = plot_from_function(p,figurenum,stretch)
    f = figure(figurenum);
    t = linspace(0,p.T,100);

    for i = 1:length(p.phi_def)
        y1(:,i) = p.phi_def{i}(t).';
    end

    if exist('stretch','var')
        y2 =  y1;
        for i = 1:100
            [y1(i,1),y1(i,2)] = stretch.old_to_new_points(y2(i,1),y2(i,2));
        end
    end
    hold on;
    plot(y1(:,1),y1(:,2),'LineWidth',3);
end