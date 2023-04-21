function subplotMN(M,N,m,n)
% Sets subplot from m,n coordinates
% ex. subplot(10,6,3,6): in a plot of 10 x 6, set plot at row 3 and column 6

if any(m>M) || any(n>N)
    error('Indexes exceed number of subplots.')
elseif any(m<0) || any(n<0)
    error('Indexes must be real integers.')
end

Idx = nan(length(m)*length(n),1);
count = 0;
for i = 1:length(m)
    for j = 1:length(n)
        count = count+1;
        Idx(count) = (m(i)-1)*N+n(j);
    end
end

subplot(M,N,Idx);
