function [E, SX, SY, SZ, SP, SM] = spinOperators(s)
    
	E = eye(2*s+1,2*s+1);
	SZ = diag([s:-1:-s]);
	SP = zeros(2*s+1,2*s+1);
	SM = zeros(2*s+1,2*s+1);
	SX = zeros(2*s+1,2*s+1);
	SY = zeros(2*s+1,2*s+1);
	for k = 1:2*s
		SP(k,k+1) = sqrt(k*(2*s+1-k));
		SM(k+1,k) = sqrt(k*(2*s+1-k));
		SX(k,k+1) = sqrt(k*(2*s+1-k))/2;
		SX(k+1,k) = sqrt(k*(2*s+1-k))/2;
		SY(k,k+1) = -1i * sqrt(k*(2*s+1-k))/2;
		SY(k+1,k) = 1i * sqrt(k*(2*s+1-k))/2;
	end

end