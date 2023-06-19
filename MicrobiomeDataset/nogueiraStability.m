function phi = nogueiraStability(sel,N)

phi = zeros(size(sel,2),1);
p = size(sel,1);

for i=1:numel(phi)
    if any(sel(:,i))
        pf = sel(:,i)/N;
        k = sum(pf);
        phi(i) = 1 - 1/(k/p * (1-k/p)) * 1/p * N/(N-1) * sum(pf.*(1-pf));
    else
        phi(i) = NaN;
    end
end