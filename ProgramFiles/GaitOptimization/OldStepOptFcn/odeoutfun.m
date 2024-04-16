function status = odeoutfun(t,y,flag,s,direction,dimension,nfparam,stretch,handles)
    status = 0;

    if strcmpi(flag,'init')
        disp('Step-optimizer starts');
        optimValues.iteration = 1;
    elseif strcmpi(flag,'done')
        disp('Step-optimizer ends');
    else
        optimValues.iteration = 2;

        global previousDisp currentDisp bestDisp;
        if isempty(previousDisp)
            previousDisp = zeros(size(currentDisp));
        end

        disp(strcat('x, y, theta: [',num2str(currentDisp.'),']'));
        ys = y(1:((nfparam-1)*dimension),end);
        ys=reshape(ys,[nfparam-1 dimension]);
        ys = [ys; 2*pi*ones(1,2)];
        
        outfun(ys,optimValues,flag,stretch,s,handles);

        % Step-optimizer stop condition
        if length(direction) == 1
            if abs(currentDisp/bestDisp) < 0.2
                status = 1;
            end
        % Steering optimizer stop condition
        else
            if norm((previousDisp(direction)-currentDisp(direction))./currentDisp(direction)) < 1e-3
                status = 1;
            end
        end

        previousDisp = currentDisp;
    end
end

