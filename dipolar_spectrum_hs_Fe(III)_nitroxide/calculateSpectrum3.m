function [spc] = calculateSpectrum3(Spin1, Spin2, Exp, Sim, Const)
%calculateSpectrum2 Calculation of the dipolar spectrum of the high-spin
%iron(III)/nitroxide spin pair via diagonalisation of the total
%Hamiltonian. The high-spin iron(III) is considered as an effective
%spin-5/2.

    % g tensor of the high-spin iron(III)
    GTM1 = diag(Spin1.gActual);
    
    % g tensor of the nitroxide
    RM2G = euler2RotationMatrix(Spin2.gFrame);
    RG2M = RM2G';
    GTM2 = RM2G * diag(Spin2.g) * RG2M;
    
    % Distance vector
    rv = spherical2Cartesian(1.0, Sim.rFrame(2), Sim.rFrame(1));
    
    % Electron spin operators
    [E1, S1X, S1Y, S1Z, S1P, S1M] = spinOperators(5/2);
    [E2, S2X, S2Y, S2Z, S2P, S2M] = spinOperators(1/2);
    
    % Direct product for electron spin operators
    SIZ = kron(S1Z,E2);
    SIP = kron(S1P,E2);
    SIM = kron(S1M,E2);
    SIX = kron(S1X,E2);
    SIY = kron(S1Y,E2);
    ISZ = kron(E1,S2Z);
    ISP = kron(E1,S2P);
    ISM = kron(E1,S2M);
    ISX = kron(E1,S2X);
    ISY = kron(E1,S2Y);
    
    % ZFS term of the Hamiltonian of the Fe(III) spin
    Dv = [-Spin1.D(1)/3 + Spin1.D(2), -Spin1.D(1)/3 - Spin1.D(2), 2*Spin1.D(1)/3];       
    DTM = diag(Dv);         
    Hzfs1 = DTM(1,1) * SIX*SIX + DTM(1,2) * SIX*SIY + DTM(1,3) * SIX*SIZ + ...
            DTM(2,1) * SIY*SIX + DTM(2,2) * SIY*SIY + DTM(2,3) * SIY*SIZ + ...
            DTM(3,1) * SIZ*SIX + DTM(2,3) * SIY*SIZ + DTM(3,3) * SIZ*SIZ;
    
    % Simulate random directions of the magnetic field
    B0v = randomDistrOnSphere(Sim.nSamples);
    
    % Calculate dipolar frequencies  
    fddv = zeros(1, Sim.nSamples);
    for k = 1:Sim.nSamples
        
        % Output progress
        if rem(k,1e5) == 0
            display(['Averaging steps ', num2str(k), '/', num2str(Sim.nSamples)])
        end
        
        % Zeeman Hamiltonian of the high-spin iron(III)
        ez1 = Const.Fez * Exp.magnField * GTM1 * B0v(:,k);
        Hez1 = (1/2) * (ez1(1) - 1i*ez1(2)) * SIP + ...
               (1/2) * (ez1(1) + 1i*ez1(2)) * SIM + ...
                ez1(3) * SIZ;       
        
        % Zeeman Hamiltonian of the nitroxide
        ez2 = Const.Fez * Exp.magnField * GTM2 * B0v(:,k);
        Hez2 = (1/2) * (ez2(1) - 1i*ez2(2)) * ISP + ...
               (1/2) * (ez2(1) + 1i*ez2(2)) * ISM + ...
                ez2(3) * ISZ;
        
        % Inter-spin distance
        r = Sim.r + Sim.rStd * randn();
        
        % Dipolar Hamiltonian
        Cdd = Const.Fdd / r^3;    
        Ddd = Cdd * ((GTM1*GTM2) - 3 * kron(GTM1*rv,(GTM2*rv)'));
        Hdd = Ddd(1,1) * SIX*ISX + Ddd(2,2) * SIY*ISY + Ddd(3,3) * SIZ*ISZ + ...
              Ddd(1,2) * SIX*ISY + Ddd(2,1) * SIY*ISX + ...
              Ddd(1,3) * SIX*ISZ + Ddd(3,1) * SIZ*ISX + ...
              Ddd(2,3) * SIY*ISZ + Ddd(3,2) * SIZ*ISY;
       
        % The total Hamiltonian
        Htot = Hzfs1 + Hez1 + Hez2 + Hdd;
        [Vtot, Wtot] = eig(Htot);
       
        % The Hamiltonians of individual spins
        H1 = Hzfs1 + Hez1;
        H2 = Hez2;
        
        % Quantum numbers
        W1 = real(Vtot' * H1 * Vtot);
        W1A = min(diag(W1)) / (-5/2);
        qn1 = diag(W1) / W1A;
        W2 = real(Vtot' * H2 * Vtot);
        W2A = min(diag(W2)) / (-1/2);

        % Order the spins states of high-spin iron(III) in accordance to its quantum numbers
        [Q, I] = sort(qn1);
        
        % Dipolar frequency
        fdd = abs( abs(Wtot(I(1),I(1)) - Wtot(I(2),I(2))) - abs(Wtot(I(3),I(3)) - Wtot(I(4),I(4))) );
        fddv(k) = fdd;
    end
    
    % Calculate the dipolar spectrum
    bins = Sim.freq;
    [spc, freq] = hist(fddv, bins);
    spc = spc + fliplr(spc);
    %spc = spc / max(spc);
    
end

