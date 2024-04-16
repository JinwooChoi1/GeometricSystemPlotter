function [eff1,eff2] = normalize_eff(eff1,eff2,anchor)

    eff1.f = (eff1.f-anchor(2,1))/(anchor(1,1)-anchor(2,1));
    eff2.f = (eff2.f-anchor(1,2))/(anchor(2,2)-anchor(1,2));

    if (isfield(eff1,'gf')&&isfield(eff2,'gf'))
        eff1.gf = eff1.gf/(anchor(1,1)-anchor(2,1));
        eff2.gf = eff2.gf/(anchor(2,2)-anchor(1,2));
    end

    if (isfield(eff1,'hf')&&isfield(eff2,'hf'))
        eff1.hf = eff1.hf/(anchor(1,1)-anchor(2,1));
        eff2.hf = eff2.hf/(anchor(2,2)-anchor(1,2));
    end
end

