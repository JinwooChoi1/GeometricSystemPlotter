function chy = chy_generator(coeff,n,dimension,invert,t0)
    %% changey/dcoeff tells us how much each point moves when a fourier series variable is changed
    % chy is a cell with as many entries as the dimension of the shape space
    % ith element of dy is a matrix where the (j,k)th entry tells us the change in the ith coordinate
    % of the kth point of the gait resulting from a unit change in the jth
    % fourier coefficient corresponding to the ith dimension of the shape space

    chy=cell(dimension,1);

    if ~exist('invert','var')
        invert = 0;
    end

    w = coeff(end,1);
    if ~exist('t0','var')
        t0 = 0;
    end

    % Create vector of time values at which to evaluate points of gait
    fourorder = length(coeff)-1;
    T = 2*pi/w;
    t = linspace(t0,T+t0,n);

    for i=1:1:dimension
        for j=1:1:n
            for k = 1:1:fourorder
                if k == 1
                    chy{i}(k,j) = 1;
                elseif mod(k,2) == 0
                    chy{i}(k,j) = cos(floor(k/2)*w*t(j));
                else
                    chy{i}(k,j) = sin(floor(k/2)*w*t(j));
                end
            end
        end
    end

    if invert == 1
        for i = 1:1:dimension
            chy{i} = flip(chy{i}.').';
        end
    end
end

