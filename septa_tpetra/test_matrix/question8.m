x(1)= 0; y(1) = -1;
h = 0.1; n = (3 -0)/h; 
f=@(y,x) (3*y*y + 2 *x*y*y); 
for i = 1:n
  x(i+1) = x(i) + h;
  y(i+1) = y(i) + h*(3 + 2 *x(i))*y(i)*y(i)
end

