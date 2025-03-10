function gait_param_to_sysplotter(y,paramfiletext)
    arguments
        y double
        paramfiletext {mustBeText}
    end

    % Load the sysplotter configuration information
    load sysplotter_config;

    %% Fourier Parameter to Gait waypoints.
    n_dim = size(y,2);
    n = 100;

    alpha_out = path_from_fourier(y,n,n_dim);

    w = y(end,1);
    t = linspace(0,2*pi/w,n+1);

    paramfiletext = strcat("opt_",paramfiletext);
    %% Save the data to a parameters file

    if (n_dim >= 2) && (n_dim <= 5)
        for i = 1:n_dim
            % The string makes "alphai = alpha_out(:,i);"
            eval(['alpha',num2str(i),'=alpha_out(:,',num2str(i),');']);
        end

        % save(fullfile(shchpath,strcat(paramfiletext,'.mat')),'alpha1','alpha2','t')
        save_string = "save(fullfile(shchpath,strcat(paramfiletext,'.mat')),";
        for i = 1:n_dim
            save_string = strcat(save_string,"'alpha",num2str(i),"',");
        end
        save_string = strcat(save_string, "'t')");
        eval(save_string);

    else
        error('Trying to make an optimal gait with an unsupported number of dimensions')
    end

    % Create the file if it doesn't already exist; future work could be
    % more fine-grained about this (making sure that any patches are
    % applied vs not overwriting any hand-edits the user made) and allowing
    % user to enter a prettier string for the display name here.

    % ['shchf_' paramfilenamebare '.m']
    filename = strcat(shchpath,'\shchf_',paramfiletext,'.m');
    if exist(filename,'file')
        delete(filename);
    end
    gait_gui_draw_make_shchf(paramfiletext,paramfiletext,n_dim);

end