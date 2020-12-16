function F = estado_estacionario(xxs)

global beta sigma delta x x_l phi phi_l phi_g lphi alpha_N alpha_T kappa omega nu d z tauu g a_N a_T xi;
c_e = xxs(1);
c = xxs(2);
l = xxs(3);
i = xxs(4);
i_N = xxs(5);
i_T = xxs(6);
k_N = xxs(7);
k_T = xxs(8);
w = xxs(9);
w_N = xxs(10);
w_T= xxs(11);
r_N = xxs(12);
r_T = xxs(13);
l_N = xxs(14);
l_T = xxs(15);
p_N = xxs(16);
y_N = xxs(17);
y_T = xxs(18);
p_X = xxs(19);
p_G = xxs(20);
b  = xxs(21);
q = xxs(22);
y = xxs(23);
lambdaa = xxs(24);
Q_N = xxs(25);
Q_T = xxs(26);
D_N = xxs(27);
r = xxs(28);
s= xxs(29);

%Sistema de ecuaciones igualadas a cero

F=[
%Ecuacion 1
c_e-((omega*(c)^((nu-1)/nu)+(1-omega)*g^((nu-1)/nu))^(nu/(nu-1)));

%Ecuacion 2
lambdaa-(omega*c^(-1/nu)*c_e^((1/nu)-1));

%Ecuacion 3
lphi*l^(sigma)-(lambdaa*(1-tauu)*w);

%Ecuacion 4
Q_N-(1+kappa*((i_N/k_N-delta))); 

%Ecuacion 5
Q_T-(1+kappa*((i_T/k_T-delta)));

%//Ecuacion_5

%Ecuacion 6
Q_N-(beta*(lambdaa/lambdaa)*((1-tauu)*r_N-kappa/2*(i_N/k_N-delta)^(2)+kappa*(i_N/k_N-delta)*(i_N/k_N)+Q_N*(1-delta)));

%Ecuacion 7
Q_T-(beta*(lambdaa/lambdaa)*((1-tauu)*r_T-kappa/2*(i_T/k_T-delta)^(2)+kappa*(i_T/k_T-delta)*(i_T/k_T)+Q_T*(1-delta)));

%Ecuacion 8
l-(((phi_l)^(-1/x_l)*(l_N)^((1+x_l)/x_l)+(1-phi_l)^(-1/x_l)*(l_T)^((1+x_l)/x_l))^(x_l/(1+x_l)));

%Ecuacion 9
l_N-(phi_l*(w_N/w)^(x_l)*l);

%Ecuacion 10
l_T-((1-phi_l)*(w_T/w)^(x_l)*l);

%Ecuacion 11
i-(i_N+i_T);

%Ecuacion 12
k_N-((1-delta)*k_N+i_N);

%Ecuacion 13
k_T-((1-delta)*k_T+i_T);

%Ecuacion 14
y_N-(a_N*k_N^(1-alpha_N)*l_N^(alpha_N));

%Ecuacion 15
alpha_N*p_N*y_N-(w_N*l_N);

%Ecuacion 16
(1-alpha_N)*p_N*y_N-(r_N*k_N);

%Ecuacion 17
(1-alpha_T)*p_X*y_T-(r_T*k_T);

%Ecuacion 18
y_T-(a_T*k_T^(1-alpha_T)*l_T^(alpha_T));

%Ecuacion 19
alpha_T*xi*s*y_T-(w_T*l_T);

%Ecuacion 20
1-((phi*(p_N)^(1-x)+(1-phi)*s^(1-x))^(1/(1-x)));

%Ecuacion 21
p_G-(phi_g*(p_N)^(1-x)+(1-phi_g)*s^(1-x))^(1/(1-x));

%Ecuacion 22
D_N-(phi*(c+i+(kappa/2)*(i_N/k_N-delta)^2*(k_N)+(kappa/2)*(i_T/k_T-delta)^(2)*k_T )+phi_g*(p_G)^(x)*g);

%Ecuacion 23
y_N-(p_N^(-x)*D_N);

%Ecuacion 24
y-(p_N*y_N+(xi*s*y_T));

%Ecuación 25
c+i+(kappa/2)*((i_N)/(k_N)-delta)^(2)*k_N+(kappa/2)*((i_T)/(k_T)-delta)^(2)*k_T+p_G*g-y-(s*(q*b-(1-d)*b));

%Ecuacion 26
q-(beta*(1-d));

%Ecuacion 27
tauu*(w*l+r_N*k_N+r_T*k_T)+q*s*b-(s*(1-d)*b+p_G*g+z);

%Ecuacion 28
r-1/q;

%Ecuacion 29
xi-p_X/s;
];