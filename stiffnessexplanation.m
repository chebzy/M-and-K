<<<<<<< HEAD
<<<<<<< HEAD
delta = 0.000001;
F = inline('y^2 - y^3','t','y');
opts = odeset('RelTol',1.e-4);
=======
delta = 0.000001;
F = inline('y^2 - y^3','t','y');
opts = odeset('RelTol',1.e-4);
>>>>>>> origin/master
=======
delta = 0.000001;
F = inline('y^2 - y^3','t','y');
opts = odeset('RelTol',1.e-4);
>>>>>>> 1465271a7a1f595c106840d24fa18e5b1f37f70c
ode45(F,[0 2/delta],delta,opts);