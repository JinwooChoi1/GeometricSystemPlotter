function options = contoptions(options)
    %%%%%%%%%%%%%
    % It is an option argument generator.
    %
    % Inputs:
    %
    % Outputs:
    %
    %%%%%%%%%%%%%

    % Default value
    arguments
        options.algorithm {mustBeMember(options.algorithm,['PredCor','PseudoArcLen'])} = 'PredCor';
        options.h (1,1) {mustBePositive} = 0.05
        options.nulltol (1,1) {mustBePositive}  = 0.2;
        options.nullHSize (1,1) {mustBeInteger,mustBeGreaterThanOrEqual(options.nullHSize,0)}  = 0;
        options.alpha (1,1) {mustBePositive}  = 0.1;
        options.breaktol (1,1) {mustBePositive}  = 0.1;
        options.corrlimit (1,1) {mustBeInteger,mustBePositive} = 100;
        options.outfcn = [];
        options.ncpfcn = [];
        options.iterlimit (1,1) {mustBeInteger,mustBePositive}  = 200;
        options.optimalitystop (1,1) {mustBePositive} = 0.1;
        options.contdir (1,1) {mustBeMember(options.contdir,[-1,1])} = 1;
    end

end