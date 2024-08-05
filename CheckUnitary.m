%% Check Unitary

function [Maxerror, S] = CheckUnitary(Filter, w)

% CheckUnitary aims to check the calculated S-parameters are unitary or not 

% w is the normalized angle frequency
% S is the S-parameters calculated by the polynomails

% By yellowbook, 2024-07-21


S(1,1,:) = polyval(Filter.Fs,w*1i)./polyval(Filter.Es,w*1i)./Filter.epsilonR;
S(1,2,:) = polyval(Filter.Ps,w*1i)./polyval(Filter.Es,w*1i)./Filter.epsilon;
S(2,1,:) = S(1,2,:);
S(2,2,:) = polyval(Filter.F22s,w*1i)./polyval(Filter.Es,w*1i)./Filter.epsilonR;
for k = 1:length(w)
    Error(k) = max(max(abs(S(:,:,k)*conj(S(:,:,k).') - eye(2))));
end
Maxerror = max(Error);  

end
