function [geff] = gTensorIron(Spin, Exp, Const)

	% Real g-tensor
    g_xx = Spin.gActual(1);
    g_yy = Spin.gActual(2);
    g_zz = Spin.gActual(3);
    gv = [g_xx, g_yy, g_zz];
    GTM = diag(gv);
    
    % Electron spin operators
    s = 5/2;
    SZ = diag([s:-1:-s]);
    SP = zeros(2*s+1, 2*s+1);
    SM = zeros(2*s+1, 2*s+1);
    SX = zeros(2*s+1, 2*s+1);
    SY = zeros(2*s+1, 2*s+1);
    for k = 1:2*s
        SP(k,k+1) = sqrt(k * (2*s+1-k));
        SM(k+1,k) = sqrt(k * (2*s+1-k));
        SX(k,k+1) = sqrt(k * (2*s+1-k))/2;
        SX(k+1,k) = sqrt(k * (2*s+1-k))/2;
        SY(k,k+1) = -k * sqrt(k * (2*s+1-k))/2;
        SY(k+1,k) = k * sqrt(k * (2*s+1-k))/2;
    end
    
    % Calculate effective g-tensor for 3 orientations along x,y,z
    orient = [0 pi/2 0; pi/2 pi/2 0; 0 0 0];
    geff = zeros([3,1]);
    for k = 1:3    
        % Rotation matrix
        RL2M = euler2RotationMatrix(orient(k,:));
        RM2L = RL2M'; 
        % Electron Zeeman term of the Hamiltonian
        GTL = RL2M * GTM * RM2L;
        GTL = (GTL + GTL')./2;
        ez = Const.Fez * GTL * [0 0 1]';
        Hez = (ez(1) - 1i*ez(2)) * SP/2 + ...
              (ez(1) + 1i*ez(2)) * SM/2 + ...
               ez(3) * SZ;
        % ZFS term of the Hamiltonian
        Dv = [-Spin.D(1)/3 + Spin.D(2), -Spin.D(1)/3 - Spin.D(2), 2*Spin.D(1)/3];       
        DTM = diag(Dv);
        DTL = RL2M * DTM * RM2L;
        DTL = (DTL + DTL')./2;         
        Hzfs = DTL(3,3) * SZ*SZ + ...
              (1/4) * (DTL(1,1) + DTL(2,2)) * SP*SM + ...
              (1/4) * (DTL(1,1) + DTL(2,2)) * SM*SP + ...
              (1/2) * (DTL(1,3) - 1i*DTL(2,3)) * (SZ*SP + SP*SZ) + ...
              (1/2) * (DTL(1,3) + 1i*DTL(2,3)) * (SZ*SM + SM*SZ) + ...
              (1/4) * (DTL(1,1) - DTL(2,2) - 1i*(DTL(1,2)+DTL(2,1))) * SP*SP + ...
              (1/4) * (DTL(1,1) - DTL(2,2) + 1i*(DTL(1,2)+DTL(2,1))) * SM*SM;
        % Calculate all resonance fields: eignenfield method
        E36 = eye(36,36);
        E6 = eye(6,6);
        A = Const.GHz2MHz * Exp.mwFreq * E36 + kron(E6, Hzfs') - kron(Hzfs, E6);
        B = -kron(E6, Hez') + kron(Hez, E6);
        b = eig(A,B);
        Nf = 0;
        for m = 1:length(b)
           if imag(b(m)) < 1e-4
               if real(b(m)) >= 0
                   if real(b(m)) < 10
                      Nf = Nf + 1;
                      Bres(Nf) = real(b(m));
                   end
               end
           end
        end
        % Find the resonance field for the lowest doublet
        B0 = 0;
        dWmin = 0.1;
        for m = 1:Nf
            Htot = Hzfs + Hez * Bres(m);
            W = eig(Htot);
            W = sort(real(W));
            dW = abs(W(1) - W(2)) - Const.GHz2MHz * Exp.mwFreq;
            if abs(dW) < dWmin
                dWmin = abs(dW);
                B0 = Bres(m);
            end           
        end
        % Calculate effective g
        geff(k) = Const.GHz2MHz * Exp.mwFreq / (Const.Fez * B0);
    end

end