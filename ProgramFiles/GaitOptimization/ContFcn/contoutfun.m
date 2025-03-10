function stop=contoutfun(y,i,s,stretch,handles,options)
    dimension = options.dimension;
    nf = options.nfparam;

    stop = false;
    optimValues.iteration = i;

    y = y(1:(nf-1)*dimension);

    y = reshape(y, [nf-1 dimension]);
    y = [y; 2*pi*ones(1,dimension)];
    outfun(y,optimValues,'',stretch,s,handles);
end