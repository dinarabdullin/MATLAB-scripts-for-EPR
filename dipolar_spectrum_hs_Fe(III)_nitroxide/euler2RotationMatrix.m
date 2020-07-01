function [RM] = euler2RotationMatrix(eulerAngles)
% Converts ZYZ Euler angles into a rotation matrix
    
    alpha = eulerAngles(1);
    beta = eulerAngles(2);
    gamma = eulerAngles(3);
    
    sa = sin(alpha);
    ca = cos(alpha);
    sb = sin(beta);
    cb = cos(beta);
    sg = sin(gamma);
    cg = cos(gamma);  
    
    RM = zeros(3,3);
    RM(1,1) = cg*cb*ca - sg*sa;
    RM(1,2) = cg*cb*sa + sg*ca;
    RM(1,3) = -cg*sb;
    RM(2,1) = -sg*cb*ca - cg*sa;
    RM(2,2) = -sg*cb*sa + cg*ca;
    RM(2,3) = sg*sb;
    RM(3,1) = sb*ca;
    RM(3,2) = sb*sa;
    RM(3,3) = cb;
    
end

