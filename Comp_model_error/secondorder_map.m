%this function maps the coordinates from reference triangle to physical
%triangle

%Input: 

%xi---(t,s) in reference triangle (quadrature point)
%V---Vertices in the physical triangle
%C---coefficient matrix of basis functions

%Output:
%Z---corresponding coorinates of 'xi' in the physical triangle




function  Z= secondorder_map(xi,V,C)
sha_fun=@(t,s,c) (c(1) + c(2)*t + c(3)*s + c(4)*t.^2 + c(5)*t.*s +c(6)*s.^2);
%compute the shape function value at (t,s) with given coefficients
n_p=size(xi,2);%number of points in the reference triangle
F_is= @(t,s,V,C) V(:,1).*sha_fun(t,s,C(:,1)) +V(:,2).*sha_fun(t,s,C(:,2)) +V(:,3).*sha_fun(t,s,C(:,3))...
              +V(:,4).*sha_fun(t,s,C(:,4)) +V(:,5).*sha_fun(t,s,C(:,5))+V(:,6).*sha_fun(t,s,C(:,6));
          %mapping to physical space

  %compute the coordinates corresponding to xi in the physical triangle        
  for i=1:n_p
      Z(:,i)=F_is(xi(1,i),xi(2,i),V,C);
  end