function y=reshape_parameter_to_twod(y,gopt)

    if iscell(y)
        for i = 1:length(y)
            y{i} = internal_function(y{i},gopt);
        end
    end

    if ismatrix(y)
        if size(y,1) == 1 || size(y,2) == 1
            y = internal_function(y,gopt);
        else
            yf = cell(size(y,1),1);
            for i = 1:size(y,1)
                yf{i} = internal_function(y(i,:).',gopt);
            end
            y = yf;
        end
    end

end

function y = internal_function(y,gopt)
    dimension = gopt.dimension;

    y = reshape(y,[], dimension);

    % If the frequency term does not exist, add it.
    if size(y,1) == gopt.nfparam - 1
        y = [y; 2*pi*ones(1,dimension)];
    end
end