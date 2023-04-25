function ZC = ZeroX(y)
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Returns Approximate Zero-Crossing Indices Of Argument Vector
zxidx = zci(y);

% the above function finds the point right before 0 crossing, not the closes one
% next step makes sure that it is the closest one

ZC=zeros(size(zxidx));
for n=1:numel(ZC)
    n0=zxidx(n);
    n1=n0+1;
    if n1<=length(y)
        if abs(y(n0))<abs(y(n1))
            ZC(n)=n0;
        else
            ZC(n)=n1;
        end
    end
end

%% remove duplicates
ZC=unique(ZC);

%% remove zero
ZC(ZC == 0) = [];

%% make sure column
ZC=ZC(:);