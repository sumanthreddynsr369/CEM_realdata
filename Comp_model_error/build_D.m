function D=build_D(n)
% %where n is the number of sine/cosine frequencies in DtN map
m=n/2;%size of NtD matrix/2
s=-[1:n]';
ss=[zeros(1,m),[1:m]]';
D=spdiags([s ss],[-m m],n,n);


