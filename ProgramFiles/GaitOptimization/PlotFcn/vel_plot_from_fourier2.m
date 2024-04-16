function vel_plot_from_fourier2(y,nfparam,figurenum)

    f = figure(figurenum);

    if (size(y,1) == 1)
        y = reshape(y.',[nfparam 2]);
    elseif(size(y,2) == 1)
        y = reshape(y,[nfparam 2]);
    end

    p = makeGait(y,0);

    %Prep figure from sysplotter
    axis tight;
    hold on;
    %colormap(cmap)

    t = linspace(0,1,50).';

    XData = p.phi_def{1}(t);
    YData = p.phi_def{2}(t);

    %Calculate distances between points on gait
    dXData = p.dphi_def{1}(t);
    dYData = p.dphi_def{2}(t);

    %Calculate length of each line segment as fraction of whole
    %gait
    percentages = sqrt(dXData.^2+dYData.^2)/...
        sum(sqrt(dXData.^2+dYData.^2));
    covariance = sqrt(cov(percentages));

    %Calculate min and max of these ratios
    minp = min(percentages);
    maxp = max(percentages);

    maxl = min(5+5e2*covariance,10);
    minl = max(5-5e2*covariance,1);

    %Set linewidths so that slower motion is a thicker line
    linewidths = interp1([minp,maxp],[maxl,minl],percentages);

    %Redraw gait with varied linewidths
    for i = 1:length(t)-1
        plot(XData(i:i+1),YData(i:i+1),'-r','LineWidth',linewidths(i));
    end

    %Draw circles under the joints between line segments to make it
    %less obvious what we did
    scatter(XData(1:length(t)),YData(1:length(t)),pi*(linewidths/2).^2,'r','filled');
end