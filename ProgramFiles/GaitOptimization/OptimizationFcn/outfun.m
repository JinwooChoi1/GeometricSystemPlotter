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

    for thisAxes = 1:numel(handles.plot_thumbnails.Children)

        ax = handles.plot_thumbnails.Children(thisAxes);

        y1 = path_from_fourier(y,n,dimension);
        hold(ax,'on');
        if stretch
            stretchnames = {'stretch','surface'};
            stretchname = stretchnames{stretch};

            [newXYZ{1},newXYZ{2},newXYZ{3}] = s.convert.(stretchname).old_to_new_points(y1(:,1),y1(:,2));
        else
            newXYZ{1} = y1(:,1);
            newXYZ{2} = y1(:,2);
            newXYZ{3} = zeros(size(y1(:,1)));
            if size(y1,2) > 2
                newXYZ{3} = y1(:,3);
            end
        end

        % The if statement below plots the gait at initial iteration
        if optimValues.iteration <= 1
            handle1=line(ax,'XData',newXYZ{1},'YData',newXYZ{2},'ZData',newXYZ{3},'color','k','linewidth',3,'UserData',{'OptimizeTracer',optimValues.iteration + 1}); %#ok<NASGU>
        end

        % The if statement below fades the gait plotted during the previous iteration
        if optimValues.iteration >= 1
            children=get(ax,'children');
            for idx = 1:numel(children)
                if iscell(children(idx).UserData) && strcmp(children(idx).UserData{1},'OptimizeTracer')
                    if children(idx).UserData{2} == 2
                        prevXYZ{1} = children(idx).XData;
                        prevXYZ{2} = children(idx).YData;
                        prevXYZ{3} = children(idx).ZData;            
                    end
                end
            end


            for idx = 1:numel(children)
                if iscell(children(idx).UserData) && strcmp(children(idx).UserData{1},'OptimizeTracer')
                    if children(idx).UserData{2} == 2
                        set(children(idx),{'Xdata','YData','ZData'},newXYZ);
                        children(idx).Color=[0 0 0];
                        children(idx).LineWidth=4;                        
                    elseif children(idx).UserData{2} == 1
                        set(children(idx),{'Xdata','YData','ZData'},prevXYZ);
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
        end
    end

    %pause(0.05)
    drawnow;
    stop=false;
end