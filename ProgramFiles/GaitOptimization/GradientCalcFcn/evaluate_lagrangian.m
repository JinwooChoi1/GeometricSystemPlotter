function lagr = evaluate_lagrangian(objct,const1,const2,A,b,p)
%%%%%%%%%%%%%
    % This function calculates the lagrangian function and its derivative
    % for step-optimizer. Some inputs are data structures containing the
    % gradient(gf) and Hessian(hf).
    %
    % Inputs:
    %
    % objct : Data structure containing the gradient and Hessian of
    %   Objective function.
    % const : Data structure containing the gradient and Hessian of
    %   Constraint function.
    %
    % Outputs:
    %
    % lagr : Data structure containing the gradient and Hessian of 
    %   the lagrangian function 
    % lambda : Lagrange multiplier is derived by the moore-penrose
    %   inverse of constraint function.
    % dlambda : The derivative of Lagrange multiplier with respect to the
    %   gait parameter.
    %%%%%%%%%%%%%

    lagr = struct();
    
    % If we consider the joint limit,
    if exist('A','var') && exist('b','var')   
        active_ineq = find(A*p-b >= 0);
        n = size(const1.gf,1);
        m = length(active_ineq);
        C = [const1.gf.'; A(active_ineq,:)];
        lambda = (C*C.')^(-1)*C*(-objct.gf);

        dlambda = zeros(1+m,n);
        for i = 1:n
            dC = [const1.hf(:,i).'; zeros(m,n)];
            CCt_inv = (C*C.')^(-1);
            dCCt_inv =  -2*CCt_inv*(dC*C.')*CCt_inv;
            dlambda(:,i) = dCCt_inv*C*(-objct.gf)+...
                CCt_inv*dC*(-objct.gf)+...
                CCt_inv*C*(-objct.hf(:,i));                          
        end


        lagr.gf = objct.gf + C.'*lambda;
        lagr.hf = objct.hf + lambda(1)*const1.hf + C.'*dlambda;
        lagr.lambda = lambda;

        if ~isempty(active_ineq)
            p_lambda = lambda;
            lambda = [];
            for i = 1:length(active_ineq)+1
                if i == 1
                    j = 1;
                    k = active_ineq(i);
                elseif i == length(active_ineq)+1
                    j = active_ineq(i-1) + 1;
                    k = size(A,1) + 1;
                else
                    j = active_ineq(i-1) + 1;
                    k = active_ineq(i);
                end
                lambda(j:k,1) = [p_lambda(i,1); zeros(k-j,1)];
            end
    
            lagr.lambda = lambda;
        end
        lagr.jointlimit = any(active_ineq);
    % Calculate the lagrange multiplier at the current fourier coefficient.
    elseif nargin == 2
        lambda = -pinv(const1.gf)*objct.gf;
        dlambda = -(1/norm(const1.gf).^2)*(objct.hf*const1.gf+const1.hf*objct.gf)...
            +(1/norm(const1.gf).^4)*(2*const1.hf*const1.gf*(const1.gf.'*objct.gf));
        lagr.gf = objct.gf + lambda*const1.gf;
        lagr.hf = objct.hf + lambda*const1.hf + const1.gf*dlambda.';
        lagr.lambda = lambda;
        lagr.jointlimit = 0;
    elseif nargin == 3
        n = size(const1.gf,1);
        A = [const1.gf.'; const2.gf.'];
        lambda = (A*A.')^-1*A*(-objct.gf);
        dlambda = zeros(2,n);
        for i = 1:n
            dA = [const1.hf(:,i).'; const2.hf(:,i).'];
            AAt_inv = (A*A.')^(-1);
            dAAt_inv =  -2*AAt_inv*(dA*A.')*AAt_inv;
            dlambda(:,i) = dAAt_inv*A*(-objct.gf)+...
                AAt_inv*dA*(-objct.gf)+...
                AAt_inv*A*(-objct.hf(:,i));
        end

        lagr.gf = objct.gf + A.'*lambda;
        lagr.hf = objct.hf + lambda(1)*const1.hf + lambda(2)*const2.hf +...
            A.'*dlambda;
        lagr.lambda = lambda;
        lagr.jointlimit = 0;
    else
        error("There might be no constraint function.")
    end
        

end



