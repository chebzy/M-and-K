L  = 6;
%cfl criteriasx
cfl =0.2;

N  = 50;
Ng = 2;
dx = L/N;
x  = 0:dx:L;
x_cells = -((Ng-1) * dx + dx/2) : dx : L + ((Ng-1) * dx + dx/2); 

m  = N + 2*Ng;

%set initial conditions
unp0 = initial_conditions(m,L, x_cells);
unp1 = unp0;

t = 0;
t_max = 2

dt = cfl * dx / max(abs(unp1));
dt = min(dt, t_max-t);


%lets Hope always
plot(x_cells,unp1);
shg
%pause(dt)  
t1 = time()
while (t < t_max)

%Initialize unp1_kg 
unp1_kg = zeros(1,N+Ng+Ng);

%apply the boundary conditions
%periodic boundary conditions
%for i = 1:Ng
%unp1(i) = unp1(N+i);
%unp1(N+Ng+i)=unp1(Ng+i);
%end
%transmissive boundary conditions
for i = 1:Ng
unp1_kg(i)     = unp1(Ng+1);
unp1_kg(N+Ng+i) = unp1(N+Ng);
end

%Apply 3rd Bc
unp1_kg(Ng+1) = unp1(N+Ng);  
  
%Define lamda,
lamda = dt/dx;

%Evaluate slope (ux) using computed cell averages, unp1 for N+1 cells
%pre-populate dummy vector to hold ux
ux = zeros(1,N+1);



%Evaluate forward and backward slopes ux
for i = Ng:N+Ng+1
  %d(i-1) = (unp1(i) - unp1(i-1))/dx;
  e(i-1) = (unp1(i+1) - unp1(i))/dx;
 end
 
for i = Ng:N+Ng+1
   d(i-1) = (unp1(i) - unp1(i-1))/dx;
end
ux = minmod(d,e);
%done



%Evaluate speed of propagation, a_p. 
u_nhalf = unp1(1:end-3) + ux(1:end-1); %11 elements(u - 1/2)
u_phalf = unp1(2:end-2) + ux(2:end);%11 elements (u + 1/2)

for i = Ng:N+1 
%u_phalf(i-1) = unp1(i-1)+ 0.5*dx *ux(i);
%u_nhalf(i-1) = unp1(i-1)- 0.5*dx *ux(i-1);
   
a_p(i-1) = max(abs(eig(u_phalf(i-1), u_nhalf(i-1))));
end
 %a_p = max(abs(eig(u_phalf, u_nhalf)));
 %a_p = max(abs(eig(u_phalf, u_nhalf)));
 a_pw = a_p(2:end); %forward half of propagation Speed (a(j+0.5))
 a_n =a_p(1:end-1); %Backward half a(j-0.5)
 
 
 %Define U_li, U_ri, ulp and urp 
 for i = Ng:N
   
U_lI(i-1) = unp1(i-1) + dx*ux(i-1)*(0.5 - lamda * a_pw(i-1));% u(j + 0.5,n)
U_rI(i-1) = unp1(i) - dx*ux(i)*(0.5 - lamda * a_pw(i-1));% u(j+1/2, n)
end
 %ux = zeros(1, N+1);

ulp = U_lI - 0.5*dt*(0.5.*(U_lI).^2);% ulp is left midpoint value of cell average 
urp =  U_rI - 0.5*dt*(0.5.*(U_rI).^2);%ulp is right midpoint value of cell average

%Flux and Flux difference
f_a = 0.5.*(ulp.^2); % flux of ulp
f_b = 0.5.*(urp.^2); % flux of urp
 %Evaluate wxp1
 flux_bac = (f_a - f_b); %flux difference ulp - urp
 flux_forw= (f_b - f_a); %flux difference urp - ulp

 %Evaluate wxp1 (w(j,n+1) and wxphalf (wj+0.5,n+0.5) at the next time step. We will need to reduce this to 51 to carter for a_n,which is 51. 
 %hence a_p becomes a_p indexed from 1 to end-1(1:end-1)(i.e a_pw) and a_n is as defined above
 deltspen = a_n - a_pw;
 deltspep = a_n + a_pw;
 jango = (lamda./(1-lamda.*deltspep));
 jango = (lamda./jango);
 kemi= (0.25*(dx -(dt.*a_pw)));
 
 

for i=Ng-1:N-1
s_timf(i) = ((1/(2*a_pw(i)))*(flux_forw(i)));
end
 
 %s_tf gives last part of equation wxphalf (w(j+0.5, n+1) s_tf is some speed * flux
 for i = Ng:N+Ng-2
 wxp1(i-1) =   unp1(i) + 0.5*dt*deltspen(i-1)*(ux(i-1)) - (jango(i-1)* flux_bac(i-1));
 %wxp1(i-1) = unp1(i) + 0.5*dt*(a_n(i-1) - a_pw(i-1))*(ux(i)) -(lamda/(1- lamda*(a_n(i-1) + a_pw(i-1)))*(flux_bac(i-1)); %new cell average based on midpoint values, at t + 1
 end
 %Evaluate wxphalf
 for i = Ng:N+Ng-2
 wxphalf(i-1) = (unp1(i) + unp1(i+1))/2 + (kemi(i-1) * (ux(i-1) - ux(i))) - s_timf(i-1); %(f_b(i) - f_a(i)));%new cell average at interface based on midpoint values, at t + 1
end  

%%%%Evaluate uxhalf (i.e ux(j+0.5)) i.e approximation of exact spatial derivatives,i.e new slope
av1= a_pw(1:end-1) - a_pw(2:end);
av1 = 1 + av1;
av2 = a_pw - a_n;
av2 = lamda * av2;
av2= 1+ av2;
av3 = a_n - a_pw;
av3 = lamda*av3;
si = 1 - av3;


%wxavf = wxp1(2:end) - wxphalf; Needed for computation of uxhalf i.e ux(j+1/2)

%Evaluate forward and backward difference for intemediate
for i = Ng: N-1 %40
wxavf(i-1) = wxp1(i) - wxphalf(i-1); %w j, n+1 - w
wxavb(i-1)= wxphalf(i-1) - wxp1(i-1); %w j+ 1/2, n+1)
end                                                                                                                                                        

p1 = wxavf./av2(1:end-1);
p2 = wxavb./ av1(1:end);

%Lets see
uxhalf = minmod(p1,p2);
uxhalf = 2/dx* uxhalf;

%Evolve to next time step % Pnew implies part of the new equation
%pnew1 = lamda*a_n(1:end-1)*wxphalf(1:end-1):

for i = Ng-1:N-1 %49
pnew1(i) = lamda*a_n(i)*wxphalf(i); %first part of equation for the final evolution of unp1, represented as unp1_kg

pnew2(i) = si(i) * wxp1(i);   %second part of equation for the final evolution of unp1

pnew3(i)= lamda * a_pw(i) *wxphalf(i); %Third part of equation for the final evolution of unp1
end

                                                                                                                                                         
for i = Ng:N-2 %47
 l_timslope1(i-1) = 0.5*dx*(lamda*a_n(i-1)).^2*uxhalf(i-1);
l_timslope2(i-1) = 0.5*dx*lamda*a_pw(i-1).^2*uxhalf(i);
end
pnew4 = l_timslope1 - l_timslope2;     %4th part of equation for the final evolution of unp1

unp2 = pnew1(1:end-2) + pnew2(1:end-2) + pnew3(1:end-2) + pnew4;

%while (t < t_max)
%initialize unp1_kg


for i = Ng+Ng:N-3 %50
unp1_kg(i) = unp2(i);
end
 %unp1 = unp1_kg;

%update time and dt
 t = t + dt;
    dt = cfl * dx / max(abs(unp1));
    dt = min(dt, t_max-t);
    lamda = dt/dx;
    plot(x_cells, unp1);
    shg
end
%unp1 = unp1_kg;
t2 = time();
s = unp1(1) + unp1(end);
s = 0.5*s;

plot(x_cells,unp1_kg, x_cells+s*t, unp0);
pause(1)
shg
disp(t2-t1)
  %compute interior fluxes
    %compute limiter phi function
 %   phi = limit(unp1,N,Ng);
  %  unp1_kg = zeros(1, N + Ng*Ng);
    
   %result = kurganov_burgers(N,unp1,t,dx,Ng,cfl,t_max,dt);
   
   %unp1_kg = result;
   
   
    
 

    
