function f = plot_jacobian_fourier(y,jacob,nfparam,dimension,figurenum)
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

    if (size(jacob,1) == 1)
        jacob = jacob.';
    end

    if size(jacob,1) == dimension*nfparam
        jacob = reshape(jacob,[nfparam dimension]);
    elseif size(jacob,1) == dimension*(nfparam - 1)
        jacob = reshape(jacob,[nfparam-1 dimension]);
        jacob = [jacob; zeros(1,dimension)];
    end

    jacob(end,:) = 2*pi*ones(1,dimension);
    jacob1=path_from_fourier(jacob,99,dimension,0);
    y1=path_from_fourier(y,99,dimension,0);
    hold on;
    quiver(y1(:,1),y1(:,2),jacob1(:,1),jacob1(:,2));
end

