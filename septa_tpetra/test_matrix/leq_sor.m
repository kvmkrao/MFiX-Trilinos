clear;

A = [3 3 4; 3 4 5; 5 2 7] 

b = [2 23 17]'

m = diag(A)

x = [ 0 0 0]'

for k=1:100
for i =1:3
    sum = 0;
    for j=1:3 
    sum = sum + A(i,j)*x(j) 
    end
    x(i) = x(i) + (b(i) - sum)/m(i); 
end 
x
end
xn = A\b
