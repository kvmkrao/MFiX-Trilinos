rho = 1664; cp = 1550; R = 8; H =12; 
Tc = 300+273; Th = 600+273; 
z =linspace(0,H,12);
T = (Th-Tc) * 0.5 * erf(z-6) + (Th + Tc)/2; 
Q = pi * R * R * rho * cp * Integral_EA2(Z, T-Tc,3); 



