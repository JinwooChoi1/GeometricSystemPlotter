function [r,t] = atan2_derivative(struct2,struct1)
    r.f = struct1.f^2 + struct2.f^2;
    r.gf = 2*(struct1.gf*struct1.f+struct2.gf*struct2.f);

    t.f = atan2(struct2.f,struct1.f);
    t.gf = (-struct2.f*struct1.gf+struct1.f*struct2.gf)/(struct1.f^2+struct2.f^2);

    if(isfield(struct1,'hf'))
        r.hf = 2*(struct1.gf*struct1.gf.'+struct1.f*struct1.hf+...
            struct2.gf*struct2.gf.'+struct2.f*struct2.hf);

        t.hf = ((struct1.f^2+struct2.f^2)*(struct2.gf*struct1.gf.'-struct2.f*struct1.hf...
            -struct1.gf*struct2.gf.'+struct2.f*struct2.hf)...
            -(2*struct1.f*struct1.gf+2*struct2.f*struct2.gf)* ...
            (-struct2.f*struct1.gf+struct1.f*struct2.gf).')...
            /(struct1.f^2+struct2.f^2)^2;
    end
end