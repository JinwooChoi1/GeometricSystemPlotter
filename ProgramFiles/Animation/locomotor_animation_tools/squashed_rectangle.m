function XY = squashed_rectangle(length,width,cap_width,points)

    % Get the X,Y points on a "squashed rectangle" as defined in the OmniGraffle
    % drawing program (a rectangle for which cap_fraction of the width is a split
    % ellipse. To visually match OmniGraffle figures, cap_fraction should be .2
    %
    % All the points returned will be on the elliptical end caps, with the
    % final points on the caps the endpoints of the straight edges. the points
    % parameter specifies how many points are included in each half-ellipse.

    %Create a parameterization for the points on the ellipse
    t = linspace(0,pi,points)';

    %Generate an ellipse whose 'a' axis is the width of the link, and whose
    % 'b' axis is twice the cap width 
    a = width/2;
    b = cap_width;

    X = -b*sin(t);
    Y = a*cos(t);

    %Offset the cap
    %X = X-(1-cap_width)*length/2;
    X = X-(length/2-cap_width);

    %Mirror the cap
    X = [X;-X;X(1)];
    Y = [Y;-Y;Y(1)];
    
    % Augment these values with a column of ones to make them amenable to
    % SE(2) transformations
    XY = [X Y ones(2*points+1,1)]';

end