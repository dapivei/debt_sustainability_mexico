clear all
clc

global beta sigma delta x x_l phi phi_l phi_g lphi alpha_N alpha_T kappa omega nu d z tauu g a_N a_T xi;

beta = 0.98;             % the discount factor


sigma = 2;               % inverse of Frisch elasti-
                                % city for labor supply   
                                
delta = 0.03;            % capital depreciation rate
                                % for capital (non-tradable
                                % and tradable sectors)
x = 0.44;                  % substitution elasticity 
                                % b/w tradables and non-
                                % tradables for ct , gt
x_l = 1;                    % substitution elasticity b/w 
                                % ltN and ltT for lt
                                
phi = 0.47;               % home bias in ct, iNt


phi_l = 0.5019;        % steady-state labor in-
                               % come share of the non-
                               % tradable sector in labor
                               % income
phi_g = .6;            
lphi = 40.6;
alpha_N = 0.6;
alpha_T = 0.55;
kappa = 1.7;
omega = .8;
nu = 0.49;
tauu = 0.227;
d = 0;
z = .03;
g = .12712;
a_N = .65;
a_T = .95;
xi = 2.1;


options = optimset('Display','iter');

xxs0 = [1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;
             1;     %20
             1;
             1;
             1;
             1;
             1;     %25
             1;
             1;
             1;
             1;
             1;     %30
             1;
             1;
             1;
             1];    %34

[xxs,Fval,exitflag] = fsolve(@estado_estacionario,xxs0,options);
exitflag

c_e = xxs(1)
c = xxs(2)
l = xxs(3)
i = xxs(4)
i_N = xxs(5)
i_T = xxs(6)
k_N = xxs(7)
k_T = xxs(8)
w = xxs(9)
w_N = xxs(10)
w_T= xxs(11)
r_N = xxs(12)
r_T = xxs(13)
l_N = xxs(14)
l_T = xxs(15)
p_N = xxs(16)
y_N = xxs(17)
y_T = xxs(18)
p_X = xxs(19)
p_G = xxs(20)
b = xxs(21)
q = xxs(22)
y = xxs(23)
lambdaa = xxs(24)
Q_N = xxs(25)
Q_T = xxs(26)
D_N = xxs(27)
r = xxs(28)
s= xxs(29)