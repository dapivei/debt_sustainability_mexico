%Fiscal limits in developing countries: A DSGE Approach. Caso: Mexico

%----------------------------------------------------------------
% 1. Variables Exogenas y Endógenas
%----------------------------------------------------------------
var 
      c_e       ${c_e}$        (long_name='consumo efectivo')
      c
      l 
      i 
      i_N        ${I^N}$
      i_T        ${I^T}$
      k_N       ${k^N}$
      k_T       ${k^T}$
      w 
      w_N      ${w^N}$
      w_T      ${w^T}$
      r_N       ${r^N}$
      r_T        ${r^T}$
      l_N 
      l_T p_N s y_N 
      y_T p_X p_G b q y lambdaa  Q_N
      Q_T D_N r tauu g a_N a_T xi DeudaPIB;

predetermined_variables 
                                        k_N k_T b;


varexo 
            v_tau    ${v_\tau}$
            v_g        ${v_g}$
            v_a_N    ${v_a^N}$
            v_a_T   ${v_a^T}$
            v_xi;

%----------------------------------------------------------------
% 2. Parámetros
%----------------------------------------------------------------

parameters 
                   beta sigma delta x x_l phi 
                   phi_l phi_g lphi alpha_N 
                   alpha_T kappa omega nu 
                   tauSS d eta gamma rho_g 
                   rho_tau rho_a rho_xi sigma_g 
                   sigma_tau sigma_a sigma_xi 
                   z xiSS a_TSS a_NSS ySS
                   gSS bSS;

beta = 0.98;                       % factor de descuento
sigma = 2;                         % inverse of Frisch elasti-
                                          % city for labor supply             
delta = 0.03;                      % capital depreciation rate
                                          % for capital (non-tradable
                                          % and tradable sectors)
x = 0.44;                            % substitution elasticity 
                                          % b/w tradables and non-
                                          % tradables for ct , gt
x_l = 1;                              % substitution elasticity b/w 
                                          % ltN and ltT for lt
phi = 0.47;                         % home bias in ct, iNt
phi_l = 0.5019;                  % steady-state labor 
                                          % income share of the 
                                          % nontradable sector 
                                          % in labor income
phi_g = 0.6;                       % home bias in gt
lphi = 40.6;                        % utility weight on leisure
alpha_N = 0.6;                  % labor income share of the
                                          % nontradable sector
alpha_T = 0.55;                 % labor income share 
                                          % of the tradable sector
kappa = 1.7;                      % investment adjustment 
                                          % cost (non-tradable and 
                                          % tradable sectors)
omega = 0.8;                     % preference weight on
                                          % ct in effective consumption
nu = 0.49;                          % elasticity of substitution b/w 
                                          % ct and gt
tauSS = 0.227;                  % income tax rate
d = 0;                                 % default decisions
eta = 0.1;                           % government spending 
                                          % response to yt−1
gamma = 0.06;                  % tau_t response to stabilize 
                                          % debt
rho_g = 0.37;                     % AR(1) coefficient in g
rho_tau = 0.77;                  % AR(1) coefficient in tau
rho_a = 0.82;                     % AR(1) coefficient in rho_a
rho_xi = 0.76;                    % AR(1) coefficient in rho_xi
sigma_g = 2.32;                % standard deviation of εg
sigma_tau = 9.11;             % standard deviation of εtau
sigma_a = 1.48;                % standard deviation of εa
sigma_xi = 2.87;               % standard deviation of εxi
z = 0.03;
xiSS = 2.1002;                  % falta verificar
a_TSS = 0.95;                 % falta verificar
a_NSS = 0.65;                % falta verificar
ySS = 0.8167;                 % falta verificar
gSS =  0.1271;                % falta verificar
bSS = .5108;                   % falta verificar

%----------------------------------------------------------------
% 3. Modelo
%----------------------------------------------------------------

%//Generar un output de las ecuaciones para LaTex
write_latex_dynamic_model;
write_latex_static_model;

model;
%Ecuaciones del problema del consumidor
%//Ecuacion_1
c_e=(omega*(c)^((nu-1)/nu)+(1-omega)*g^((nu-1)/nu))^(nu/(nu-1));
%//Ecuacion_2
lambdaa=omega*c^(-1/nu)*c_e^((1/nu)-1);
%//Ecuacion_3
lphi*l^(sigma)=lambdaa*(1-tauu)*w;
%//Ecuacion_4
Q_N=1+kappa*((i_N/k_N-delta));
%//Ecuacion_5
Q_T=1+kappa*((i_T/k_T-delta));
%//Ecuacion_6
Q_N=beta*(lambdaa(+1)/lambdaa)*((1-tauu(+1))*r_N(+1)-kappa/2*(i_N(+1)/k_N(+1)-delta)^(2)+kappa*(i_N(+1)/k_N(+1)-delta)*(i_N(+1)/k_N(+1))+Q_N(+1)*(1-delta));
%//Ecuacion_7
Q_T=beta*(lambdaa(+1)/lambdaa)*((1-tauu(+1))*r_T(+1)-kappa/2*(i_T(+1)/k_T(+1)-delta)^(2)+kappa*(i_T(+1)/k_T(+1)-delta)*(i_T(+1)/k_T(+1))+Q_T(+1)*(1-delta));
%//Ecuacion_8
l=((phi_l)^(-1/x_l)*(l_N)^((1+x_l)/x_l)+(1-phi_l)^(-1/x_l)*(l_T)^((1+x_l)/x_l))^(x_l/(1+x_l));
%//Ecuacion_9
l_N=phi_l*(w_N/w)^(x_l)*l;
%//Ecuacion_10
l_T=(1-phi_l)*(w_T/w)^(x_l)*l ;
%//Ecuacion_11
i=i_N+i_T;
%//Ecuacion_12
k_N(+1)=(1-delta)*k_N+i_N;
%//Ecuacion_13
k_T(+1)=(1-delta)*k_T+i_T;

%Ecuaciones del problema del productor
%//Ecuacion_14
y_N=a_N*k_N^(1-alpha_N)*l_N^(alpha_N);
%//Ecuacion_15
alpha_N*p_N*y_N=w_N*l_N;
%//Ecuacion_16
(1-alpha_N)*p_N*y_N=r_N*k_N;
%//Ecuacion_17
(1-alpha_T)*p_X*y_T=r_T*k_T;
%//Ecuacion_18
y_T=a_T*k_T^(1-alpha_T)*l_T^(alpha_T);
%//Ecuacion_19
alpha_T*xi*s*y_T=w_T*l_T;
%//Ecuacion_20
1=(phi*(p_N)^(1-x)+(1-phi)*(s)^(1-x))^(1/(1-x));


%Ecuaciones del problema del Gobierno
%//Ecuacion_21
p_G=(phi_g*(p_N)^(1-x)+(1-phi_g)*(s)^(1-x))^(1/(1-x));
%//Ecuacion_22
D_N=phi*(c + i + (kappa/2)*(i_N/k_N - delta)^2*(k_N)+(kappa/2)*(i_T/k_T-delta)^(2)*k_T) + phi_g*(p_G)^(x)*g;
%//Ecuacion_23
y_N=p_N^(-x)*D_N;
%//Ecuacion_24
y=p_N*y_N + xi*s*y_T;
%//Ecuacion_25
c+i+(kappa/2)*((i_N)/(k_N)-delta)^(2)*k_N + (kappa/2)*((i_T)/(k_T)-delta)^(2)*k_T + p_G*g - y = s*(q*b(+1)-(1-d)*b);
%//Ecuacion_26
q=beta*(1-d);
%//Ecuacion_27
tauu*(w*l+r_N*k_N+r_T*k_T)+q*s*b(+1)=s*(1-d)*b+p_G*g+z;
%Procesos estocásticos
%//Ecuacion_28
log(tauu/tauSS)= (gamma/(1-rho_tau))*log((1-d)*b/bSS)+ v_tau;
%//Ecuacion_29
log(g/gSS)= rho_g*log(g(-1)/gSS)+eta*log(y(-1)/ySS)+ v_g;
%//Ecuacion_30
log(a_N/a_NSS)= rho_a*log(a_N(-1)/a_NSS)+ v_a_N;
%//Ecuacion_31
log(a_T/a_TSS)= rho_a*log(a_T(-1)/a_TSS)+ v_a_T;
%//Ecuacion_32
log((xi)/xiSS)= rho_xi*log(xi(-1)/xiSS)+ v_xi;
%Otros
%//Ecuacion_33
r=1/q;
%//Ecuacion_34
xi=p_X/s;
%//Ecuacion_34
DeudaPIB=b/y;
end;


%----------------------------------------------------------------
% 4. Estado estacionario
%----------------------------------------------------------------

initval;
c_e = 0.3139;
c = 0.505;
l = 0.2062;
i = 0.1563;
i_N = 0.102;
i_T = 0.0544;
k_N = 3.3993;
k_T = 1.8118;
w = 2.3128;
w_N = 2.7259;
w_T = 1.8032;
r_N = 0.0652;
r_T = 0.0652;
l_N = 0.122;
l_T = 0.0801;
p_N = 1.8467;
y_N = 0.3001;
y_T = 0.3096;
p_X = 0.848;
p_G = 1.1866;
b = 0.5108;
q = 0.98;
y = 0.8167;
lambdaa = 0.9656;
Q_N = 1;
Q_T = 1;
D_N = 0.3931;
r = 1.0204;
s = 0.4463;
g = 0.12712;
tauu = 0.227;
a_N = .65;
a_T = .95;
xi= 1.90;
DeudaPIB = .63;
end;

model_diagnostics;

steady;                              % verificar que se trata 
                                         %  de un estado estacionario
check;                               %  verificar las condiciones 
                                         %  de Blanchar-Kahn

%----------------------------------------------------------------
% 5. Choques
%----------------------------------------------------------------
shocks;
var v_tau; stderr 9.11;
var v_g; stderr 2.32;
var v_a_N; stderr 1.48;
var v_a_T; stderr 1.48;
var v_xi; stderr 2.87;
end;
colormap winter;
stoch_simul(periods=500, order=1);

%----------------------------------------------------------------
% 6. Choques simultaneos
%----------------------------------------------------------------

varobs y;
shock_decomposition(parameter_set=calibration,datafile=data) y c_e g DeudaPIB;
