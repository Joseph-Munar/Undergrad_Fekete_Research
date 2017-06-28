function DBr = DBrClosed(N,K)
%Maps local derivatives to global derivatives
% constucted in closed form (for numerical stability)
n = N*K/2;
a = min([N K]);
if N == K
    a = N-1;
end
DBr = zeros(N+K-1,N+K);
for i = 1:length(DBr(:,1))-a
    if i < a
        DBr(i,i:i+1) = [-n/i, n/i];
        DBr(end-i+1, end-i:end-i+1) = [-n/i, n/i];
    else
        DBr(i,i:i+1) = [-n/a, n/a];
        DBr(end-i+1, end-i:end-i+1) = [-n/a, n/a];
    end
end
return