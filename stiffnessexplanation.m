delta = 0.01;
F = inline('y^2 - y^3','t','y');
opts = odeset('RelTol',1.e-4);
[t,y] = ode45(F,[0, 100],delta,opts);
plot(t,y)