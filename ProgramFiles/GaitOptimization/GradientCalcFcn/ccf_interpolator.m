function [ccf,dccf] = ccf_interpolator(y,s,n,dimension,direction)
% Preallocating memory for variables which we will need in further
% calculation 
interpstateccf=cell(1,dimension); % Variable which will store the ccf function grid used for interpolation
ccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store ccf function at each point
dccf=zeros(n,dimension*(dimension-1)/2); % Variable which will store the gradient of ccf function at each point

% Interpolation to calculate all the variables needed for gradient
% calculation
y_midpoint = zeros(size(y));
y_midpoint(1,:) = (y(n,:)+y(1,:)+y(2,:))/3;
for i = 2:n-1
    y_midpoint(i,:) = (y(i-1,:)+y(i,:)+y(i+1,:))/3;
end
y_midpoint(n,:) = (y(n-1,:)+y(n,:)+y(1,:))/3;

y_for_interp = mat2cell(y_midpoint,size(y_midpoint,1),ones(1,size(y_midpoint,2)));

for j=1:1:dimension
    interpstateccf{j}=s.grid.eval{j,1};
end

for j=1:dimension*(dimension-1)/2
    ccf(:,j)=interpn(interpstateccf{:},s.DA_optimized{direction,j},y_for_interp{:},'spline');
end


% dccf calculation
if nargout == 2
    if dimension ==2
        dDA_optimized = cell(2,1);
        dr = (s.grid_range(2)-s.grid_range(1))/size(s.DA_optimized{direction},1);
        [dDA_optimized{1}, dDA_optimized{2}] = gradient(s.DA_optimized{direction});
        dDA_optimized{1} = dDA_optimized{1}/dr;
        dDA_optimized{2} = dDA_optimized{2}/dr;
    
        for j=1:dimension
            dccf(:,j)=interpn(interpstateccf{:},dDA_optimized{j},y_for_interp{:},'spline');
        end
    end
end


% ddccf calculation
if nargout == 3
    if dimension ==2
        ddDA_optimized = cell(2,1);
        gradient(s.DA_optimized{direction});
    
        for j=1:dimension
            dccf(:,j)=interpn(interpstateccf{:},dDA_optimized{j},y_for_interp{:},'spline');
        end
    end
end

end