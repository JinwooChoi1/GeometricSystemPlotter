function [X,dX] = fischerburmfunc(mui,gi,dgi)
    % It gives the result of the Fischer-Burmeister function
    % and the Jacobian of it. This function has two argument. 
    %
    % For the complementary problem,
    % mu : The Lagrange multiplier for the inequality
    % g : The inequality function.
    % dg : The gradient of the inequality function (nx1 vector).

    X = sqrt(mui^2+gi^2) - mui + gi;
    if exist('dgi','var')
        dX{1} = mui/sqrt(mui^2+gi^2) - 1;
        dX{2} = dgi*(gi/sqrt(mui^2+gi^2) + 1);
    else
        dX = [];
    end
end