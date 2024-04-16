function newstruct = fraction_derivative(struct1,struct2)
    % struct1.f, struct2.f: a scalar value.
    % struct1.gf, struct2.gf: a gradient of struct.
    % Calculate d^2(struct1.f/struct2.f).
    
    newstruct = struct();

    if(isfield(struct1,'hf'))
        newstruct.hf = (struct1.hf/struct2.f)-...
            (2/struct2.f^2)*(struct1.gf)*(struct2.gf.')-...
            (struct1.f/struct2.f^2)*(struct2.hf)-...
            2*(struct1.f/struct2.f^3)*(struct2.gf)*(struct2.gf.');
    end

    newstruct.gf = struct1.gf/struct2.f-struct1.f/struct2.f^2*struct2.gf;

    newstruct.f = struct1.f/struct2.f;
end