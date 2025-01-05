%this function helps to find coefficients of shape function
%Ac=I
%c--coefficients of shape function
%this function gives A

function row=p_eval(t,s)
f1=@(t,s) 1;
f2=@(t,s) t;
f3=@(t,s) s;
f4=@(t,s) t^2;
f5=@(t,s) t*s;
f6=@(t,s) s^2;
row=[f1(t,s),f2(t,s),f3(t,s),f4(t,s),f5(t,s),f6(t,s)];

