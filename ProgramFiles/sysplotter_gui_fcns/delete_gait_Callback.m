% --- Executes on button press in delete_gait.
function delete_gait_Callback(hObject, eventdata, handles)
% hObject    handle to delete_gait (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    

    % Update the waitbar to indicate that the process has started
    waitbar2a(0,handles.progresspanel,'Find the corresponding data');
    
    shch_index = get(handles.shapechangemenu,'Value');
    shch_names = get(handles.shapechangemenu,'UserData');
    
    current_shch = shch_names{shch_index};
    current_shch = replace(current_shch,'shchf_','');

    % Get the setup configuration file
    configfile = './sysplotter_config';
    load(configfile,'datapath','shchpath');

    % Get the filelist in datapath and shchpath
    shchfilelist = what(shchpath);
    datafilelist = what(datapath);

    % Find the file which match with the selected shape chages
    deleteshchmlist = find(contains(shchfilelist.m,current_shch));
    deleteshchmatlist = find(contains(shchfilelist.mat,current_shch));
    deletedatamlist = find(contains(datafilelist.m,current_shch));
    deletedatamatlist = find(contains(datafilelist.mat,current_shch));

    %Show half progress bar
    waitbar2a(0.5,handles.progresspanel,'Deleting');

    % If the matching file is found, delete the file.
    if ~isempty(deleteshchmlist)
        for i = deleteshchmlist.'
            deletefullfile = fullfile(shchpath, shchfilelist.m{i});
            delete(deletefullfile);
        end
    end

    if ~isempty(deleteshchmatlist)
        for i = deleteshchmatlist.'
            deletefullfile = fullfile(shchpath, shchfilelist.mat{i});
            delete(deletefullfile);
        end
    end

    if ~isempty(deletedatamlist)
        for i = deletedatamlist.'
            deletefullfile = fullfile(datapath, datafilelist.m{i});
            delete(deletefullfile);
        end
    end


    if ~isempty(deletedatamatlist)
        for i = deletedatamatlist.'
            deletefullfile = fullfile(datapath, datafilelist.mat{i});
            delete(deletefullfile);
        end
    end
    
    % Refresh the shape change selection.
    refresh_gui_Callback(hObject,eventdata,handles);

    % Show full progress bar
    waitbar2a(1,handles.progresspanel,'Finished deleting');
end