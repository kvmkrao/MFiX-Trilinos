A = pi* 2.5 * 2.5; dt = 20; R = 5;  
f = @(H) -0.55 *A*sqrt(2*9.81*H);
H(1) = 2; t(1) = 0; n = 500/dt; v(1) = 0 
for i = 1:n
    H(i+1) = H(i) + f(H(i))*dt;
    t(i+1) = t(i) + dt; 
    v(i+1) = pi * R * R *(3 * 3.6/2 - H(i))/3
end



