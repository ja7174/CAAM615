%% CAAM 615 Final Project
%% Author: Javeria Ahmed

function H_sub = HinesSub(N)
H_sub = diag(-2*ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
H_sub(1,1)=-1;
end