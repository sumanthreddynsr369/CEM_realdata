
%This is used in the exact propagation of NtD matrix  
function D=ydiag(t,n)
  %n-number of sine/cosine frequencies
  %t- p_0/p_1 - at specific ring of disc
  
  for j=1:n
      s(j)=t^(2*j);
  end
  D=diag(s);
end