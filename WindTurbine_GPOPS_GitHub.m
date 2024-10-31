% function [t_gpops,wind_profile_gpops,theta_gpops,Pin_gpops,DthetaDt_gpops] = GreenTech_GPOPS(wind_profile_type,mu_gpops,std_gpops,step_value_gpops,final_control,final_states)
%-------------------------- Wind-Turbine Optimal control Problem --------------------------%
clc;
clear all
% format long

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%

t0 = 0;
tf = 70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter the parameters below for the wind profile
b=11;
m=10;
t_on=0.2;
t_off=0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REF=[b m t_on t_off]; % reference parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finding steady state to use as initial value for the simulations
v0 = b;
u0 = 0;
options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
Yo = 1*ones(10,1);
g= @(y) steady_states(y,u0,v0); 
[y_sol,fval] = fsolve(g,Yo,options_steady_states); % finding equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_states = y_sol';
ymin = [0.1,0.1,-2,-3,0.7,0.7,0.9,0.5,0.5,0.7];
ymax = [0.4,0.3,2,3,1.5,1.3,1.1,1.45,1.45,1.45];
x0 = y_states;
xf = x0;
xMin = ymin;
xMax = ymax;
uMin = 0;
uMax = 26;
betamin = -0.3;
betamax = 0.5;
%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;
bounds.phase.initialstate.lower = x0; 
bounds.phase.initialstate.upper = x0;
bounds.phase.state.lower = xMin;
bounds.phase.state.upper = xMax;
bounds.phase.finalstate.lower = xMin; 
bounds.phase.finalstate.upper = xMax;
bounds.phase.control.lower = uMin; 
bounds.phase.control.upper = uMax;
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 100000;
bounds.parameter.lower = betamin;
bounds.parameter.upper = betamax;
%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
betaguess = 0;
guess.phase.time  = [t0; tf]; 
guess.phase.state = [[x0(1);xf(1)],[x0(2);xf(2)],[x0(3);xf(3)],[x0(4);xf(4)],[x0(5);xf(5)],[x0(6);xf(6)],[x0(7);xf(7)],[x0(8);xf(8)],[x0(9);xf(9)],...
                     [x0(10);xf(10)]];
guess.phase.control = [1.19520455299623; 1.19520455299623];
guess.phase.control = [u0; u0];
guess.phase.integral = 0;
guess.parameter = [betaguess];
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%

mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-6;
mesh.maxiterations   = 4;
mesh.colpointsmin    = 6;
mesh.colpointsmax    = 6;
mesh.phase.colpoints = 4*ones(1,10);
mesh.phase.fraction  = 0.1*ones(1,10);


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mu_gpops = 20;
std_gpops = 0.5;
auxdata.v0 = b;
auxdata.mu_gpops = mu_gpops;
auxdata.std_gpops = std_gpops;
auxdata.ramp_slope = (12.5-11.5)/(21-19);
auxdata.wind_profile_type = 'Gaussian';
% auxdata.objfun = 'abs';
auxdata.objfun = 'min';
%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'Wind-Turbine-Problem';
setup.functions.continuous           = @WindTurbineContinuous_GPOPS;
setup.functions.endpoint             = @WindTurbineEndpoint_GPOPS;
setup.displaylevel                   = 2;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'ipopt';
setup.nlp.snoptoptions.tolerance     = 1e-9;
setup.nlp.snoptoptions.maxiterations = 20000;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-10;
setup.derivatives.supplier           = 'sparseCD';

setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-bounds';

%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc

% Collecting results
result_states = output.result.solution.phase.state;
result_control = output.result.solution.phase.control;
result_states_control = [output.result.solution.phase.state,output.result.solution.phase.control];
time_opt = output.result.solution.phase.time;
theta_opt = result_control;
%% 
R1=0.02;
X1=0.0243;
Xtr=0.00557;
E=1.0164;
R=R1;
X=X1+Xtr;
Xeq=0.8;

A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
B = -(2.*result_states(:,10).*R + 2*X.*result_states(:,9)./Xeq + 2*(R^2+X^2).*result_states(:,9)./Xeq);
C = (R^2+X^2)/Xeq + (R^2+X^2).*result_states(:,10).^2 - E^2;
V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);

w0=1;
T_pl=0.3;
k_pp=150;
k_ip=25;
k_pc=3;
k_ic=30;
D_tg=1.5;
k_tg=1.11;
P_stl=1;

Tpwr=0.05;
Kqi=0.1;
Kvi=40;
Tpc=0.05;
Kitrq=0.6;
Kptrq=3;
H=4.33;
Hg=0.62;
Ktg=1.11;
Dtg=1.5;
wb=125.66;
k_b=56.6;
Kb=k_b;
Cmech=0.00159; %1/2*rho*A_r

Pelec_opt = result_states(:,10).*V;
if strcmp(auxdata.wind_profile_type,'Gaussian')
    rand_profile = (1/(std_gpops*sqrt(2*pi)))*exp(-1/2*((time_opt-mu_gpops)/std_gpops).^2);
    wind_profile = b + rand_profile;
elseif strcmp(auxdata.wind_profile_type,'ramp')
    wind_profile = b;
    wind_step = wind_profile;
    ramp_slope = auxdata.ramp_slope;
    for tt = 2:length(time_opt)
        if (time_opt(tt) < 19) || (time_opt(tt) > 21)
            ramp_profile = 0;
        else
            ramp_profile = ramp_slope*(time_opt(tt)-time_opt(tt-1));
            wind_step = wind_step + ramp_profile;
        end
     wind_profile = [wind_profile; wind_step];
    end
end
P_mech=((Cmech*wind_profile.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(result_states(:,2)+w0)./wind_profile)-1.2406e-002*(k_b.*(result_states(:,2)+w0)./wind_profile).^2-1.3365e-004*(k_b.*(result_states(:,2)+w0)./wind_profile).^3+1.1524e-005*(k_b.*(result_states(:,2)+w0)./wind_profile).^4 ...
                         -6.7606e-002.*result_control +6.0405e-002.*result_control.*(k_b.*(result_states(:,2)+w0)./wind_profile)     -1.3934e-002.*result_control.*(k_b.*(result_states(:,2)+w0)./wind_profile).^2    +1.0683e-003.*result_control.*(k_b.*(result_states(:,2)+w0)./wind_profile).^3    -2.3895e-005.*result_control.*(k_b.*(result_states(:,2)+w0)./wind_profile).^4 ...
                         +1.5727e-002.*result_control.^2  -1.0996e-002.*result_control.^2.*(k_b.*(result_states(:,2)+w0)./wind_profile)   +2.1495e-003.*result_control.^2.*(k_b.*(result_states(:,2)+w0)./wind_profile).^2  -1.4855e-004.*result_control.^2.*(k_b.*(result_states(:,2)+w0)./wind_profile).^3  +2.7937e-006.*result_control.^2.*(k_b.*(result_states(:,2)+w0)./wind_profile).^4 ...
                         -8.6018e-004.*result_control.^3  +5.7051e-004.*result_control.^3.*(k_b.*(result_states(:,2)+w0)./wind_profile)   -1.0479e-004.*result_control.^3.*(k_b.*(result_states(:,2)+w0)./wind_profile).^2  +5.9924e-006.*result_control.^3.*(k_b*(result_states(:,2)+w0)./wind_profile).^3  -8.9194e-008.*result_control.^3.*(k_b.*(result_states(:,2)+w0)./wind_profile).^4 ...
                         +1.4787e-005.*result_control.^4  -9.4839e-006.*result_control.^4.*(k_b.*(result_states(:,2)+w0)./wind_profile)   +1.6167e-006.*result_control.^4.*(k_b.*(result_states(:,2)+w0)./wind_profile).^2  -7.1535e-008.*result_control.^4.*(k_b*(result_states(:,2)+w0)./wind_profile).^3  +4.9686e-010.*result_control.^4.*(k_b.*(result_states(:,2)+w0)./wind_profile).^4))));
% Plot results
fpath = 'C:\Users\hesha\Desktop\OneDrive\OneDrive - University of Cincinnati\MATLAB\WindTurbineGPOPS2\RL Project\figs';
configuration_text = ['u0=',num2str(u0),' vw=',num2str(b),' std=',num2str(std_gpops), ' GPOPS'];

close all
set(0,'DefaultFigureWindowStyle','docked')
font_size = 20;
Pmech_fig_GPOPS = figure('Name','Pmech','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
theta_fig_GPOPS = figure('Name','u(t)','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
v_fig_GPOPS = figure('Name','v_w','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
linecolor = 'k';

figure(Pmech_fig_GPOPS)
hold on
plot(time_opt,P_mech,linecolor)
hold on
xlabel('t')
ylabel('Pmech')
title([configuration_text]);
fname = ['Pmech, ',configuration_text,'.png'];
xlim([0 60]);
% ylim([0.8 1.2]);
saveas(gca, fullfile(fpath, fname));

figure(theta_fig_GPOPS)
hold on
plot(time_opt,theta_opt,linecolor)
hold on
xlabel('t')
ylabel('u(t)')
title([configuration_text])
fname = ['u(t), ',configuration_text,'.png'];
xlim([0 60]);
ylim([0 4]);
saveas(gca, fullfile(fpath, fname));

figure(v_fig_GPOPS)
plot(time_opt,wind_profile,linecolor)
xlabel('t')
ylabel('v_w')
title([configuration_text])
fname = ['v_w, ',configuration_text,'.png'];
xlim([0 60]);
saveas(gca, fullfile(fpath, fname));


% end

function f = steady_states(y,u0,v0)
w0=1;
T_pl=0.3;
k_pp=150;
k_ip=25;
k_pc=3;
k_ic=30;
D_tg=1.5;
k_tg=1.11;
P_stl=1;
%%%%%%%%%%%%


Tpwr=0.05;
Kqi=0.1;
Kvi=40;
Xeq=0.8;
Tpc=0.05;
Kitrq=0.6;
Kptrq=3;
H=4.33;
Hg=0.62;
Ktg=1.11;
Dtg=1.5;
wb=125.66;
k_b=56.6;
Cmech=0.00159; %1/2*rho*A_r


PFE_ref=pi/20;

v_w=v0;

u=u0;
%******************************
% syms V_term Q Eq Ip
% P=V_term*Ip;
% Q=V_term*(Eq-V_term)/Xeq;
R1=0.02;
X1=0.0243;
Xtr=0.00557;
E=1.0164;
R=R1;
X=X1+Xtr;

A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
B = -(2.*y(10).*R + 2*X.*y(9)./Xeq + 2*(R^2+X^2).*y(9)./Xeq);
C = (R^2+X^2)/Xeq + (R^2+X^2).*y(10).^2 - E^2;
V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);

% V=y(13);

f=[((-V*y(10))/(y(1)+1)-y(3)*Ktg-Dtg*(y(1)-y(2)))/(2*Hg);    
((Cmech*v_w^3*((-4.1909e-001+2.1808e-001*(k_b*(y(2)+w0)/v_w)-1.2406e-002*(k_b*(y(2)+w0)/v_w)^2-1.3365e-004*(k_b*(y(2)+w0)/v_w)^3+1.1524e-005*(k_b*(y(2)+w0)/v_w)^4 ...
                              -6.7606e-002*u    +6.0405e-002*u*(k_b*(y(2)+w0)/v_w)     -1.3934e-002*u*(k_b*(y(2)+w0)/v_w)^2    +1.0683e-003*u*(k_b*(y(2)+w0)/v_w)^3    -2.3895e-005*u*(k_b*(y(2)+w0)/v_w)^4 ...
                              +1.5727e-002*u^2  -1.0996e-002*u^2*(k_b*(y(2)+w0)/v_w)   +2.1495e-003*u^2*(k_b*(y(2)+w0)/v_w)^2  -1.4855e-004*u^2*(k_b*(y(2)+w0)/v_w)^3  +2.7937e-006*u^2*(k_b*(y(2)+w0)/v_w)^4 ...
                              -8.6018e-004*u^3  +5.7051e-004*u^3*(k_b*(y(2)+w0)/v_w)   -1.0479e-004*u^3*(k_b*(y(2)+w0)/v_w)^2  +5.9924e-006*u^3*(k_b*(y(2)+w0)/v_w)^3  -8.9194e-008*u^3*(k_b*(y(2)+w0)/v_w)^4 ...
                              +1.4787e-005*u^4  -9.4839e-006*u^4*(k_b*(y(2)+w0)/v_w)   +1.6167e-006*u^4*(k_b*(y(2)+w0)/v_w)^2  -7.1535e-008*u^4*(k_b*(y(2)+w0)/v_w)^3  +4.9686e-010*u^4*(k_b*(y(2)+w0)/v_w)^4))...
          /(y(2)+w0))+D_tg*(y(1)-y(2))+k_tg*y(3))/(2*H);
wb*(y(1)-y(2));
y(1)+1-1.2;
% y(7)-P_stl;
% (k_pp*(y(1)+w0-1.2)+k_ip*y(4)+k_pc*(y(7)-P_stl)+k_ic*y(5)-y(6))/T_pl;
((y(1)+1)*(Kptrq*(y(1)+1-1.2) + Kitrq*y(4))-y(5))/Tpc;
(V*y(8)-y(6))/Tpwr;
(tan(PFE_ref)*y(6)-(V)*(y(9)-V)/Xeq)*Kqi;
(y(7)-V)*Kvi;
(y(8)-y(9))/0.02;
((y(5)/V)-y(10))/0.02];
% Equation of y13
% (y(13))^4-(2*(y(13)*y(12)*R+y(13)*(y(11)-y(13))/Xeq*X)+E^2)*(y(13))^2+(R^2+X^2)*((y(13)*y(12))^2+(y(13)*(y(11)-y(13))/Xeq)^2)];


end

%---------------------------------%
% BEGIN: WindTurbineContinuous_GPOPS.m %
%---------------------------------%
function phaseout = WindTurbineContinuous_GPOPS(input)

% format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter the parameters below for the wind profile
b=input.auxdata.v0;
m=10;
t_on=0.2;
t_off=0.3;


t = input.phase.time;
p = input.phase.parameter;
beta = p(:,1);

if isfield(input,'auxdata')
    mu = input.auxdata.mu_gpops;
    std = input.auxdata.std_gpops;
    wind_profile_type = input.auxdata.wind_profile_type;
    ramp_slope = input.auxdata.ramp_slope;
    if strcmp(wind_profile_type,'Gaussian')
        rand_profile = (1/(std*sqrt(2*pi)))*exp(-1/2*((t-mu)/std).^2);
        wind_profile = b + rand_profile;      
    elseif strcmp(wind_profile_type,'ramp')
        wind_profile = b;
        wind_step = wind_profile;
        for tt = 2:length(t)
            if (t(tt) < 19) || (t(tt) > 21)
                ramp_profile = 0;
            else
                ramp_profile = ramp_slope*(t(tt)-t(tt-1));
                wind_step = wind_step + ramp_profile;
            end
         wind_profile = [wind_profile; wind_step];
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REF=[b m t_on t_off]; % reference parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0=1;
T_pl=0.3;
k_pp=150;
k_ip=25;
k_pc=3;
k_ic=30;
D_tg=1.5;
k_tg=1.11;
P_stl=1;
%%%%%%%%%%%%


Tpwr=0.05;
Kqi=0.1;
Kvi=40;
Tpc=0.05;
Kitrq=0.6;
Kptrq=3;
H=4.33;
Hg=0.62;
Ktg=1.11;
Dtg=1.5;
wb=125.66;
k_b=56.6;
Kb=k_b;
Cmech=0.00159; %1/2*rho*A_r


b=REF(1);
m=REF(2);
t_on=REF(3);
t_off=REF(4);

V_w = wind_profile;
% v_w=V_w;
PFE_ref=pi/20;
 
%******************************
% syms V_term Q Eq Ip
% P=V_term*Ip;
% Q=V_term*(Eq-V_term)/Xeq;

%******************************
y = input.phase.state;
% y(:,5)=beta+1;

v_w = V_w;
u = input.phase.control;

R1=0.02;
X1=0.0243;
Xtr=0.00557;
E=1.0164;
R=R1;
X=X1+Xtr;
Xeq=0.8;

A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
B = -(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq);
C = (R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2;
V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);


dy_dt(:,1)=(-((-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq))).*y(:,10))./(y(:,1)+1)-y(:,3).*Ktg-Dtg.*(y(:,1)-y(:,2))./(2*Hg);    
dy_dt(:,2)=((Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(y(:,2)+w0)./v_w)-1.2406e-002*(k_b.*(y(:,2)+w0)./v_w).^2-1.3365e-004*(k_b.*(y(:,2)+w0)./v_w).^3+1.1524e-005*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                              -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(y(:,2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(y(:,2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(y(:,2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                              +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(y(:,2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                              -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(y(:,2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(y(:,2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(y(:,2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                              +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(y(:,2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(y(:,2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(y(:,2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(y(:,2)+w0)./v_w).^4))...
          ./(y(:,2)+w0))+D_tg.*(y(:,1)-y(:,2))+k_tg*y(:,3))./(2*H);
dy_dt(:,3)=wb.*(y(:,1)-y(:,2));
dy_dt(:,4)=y(:,1)+1-1.2;
dy_dt(:,5)=((y(:,1)+1).*(Kptrq.*(y(:,1)+1-1.2) + Kitrq.*y(:,4))-y(:,5))./Tpc;
dy_dt(:,6)=(((-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq))).*y(:,8)-y(:,6))./Tpwr;
dy_dt(:,7)=(tan(PFE_ref).*y(:,6)-((((-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq))).*(y(:,9)-(-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)))./Xeq)))*Kqi;
dy_dt(:,8)=(y(:,7)-(-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq))).*Kvi;
dy_dt(:,9)=(y(:,8)-y(:,9))./0.02;
dy_dt(:,10)=((y(:,5)./((-1.*(-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq))+sqrt((-(2.*y(:,10).*R + 2*X.*y(:,9)./Xeq + 2*(R^2+X^2).*y(:,9)./Xeq)).^2-4*(1 + 2*X/Xeq + (R^2+X^2)/Xeq)*((R^2+X^2)/Xeq + (R^2+X^2).*y(:,10).^2 - E^2)))./(2*(1 + 2*X/Xeq + (R^2+X^2)/Xeq))))-y(:,10))/0.02;
 

P_stl = 1;
P_inp = y(:,10).*V;
P_mech=((Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(y(:,2)+w0)./v_w)-1.2406e-002*(k_b.*(y(:,2)+w0)./v_w).^2-1.3365e-004*(k_b.*(y(:,2)+w0)./v_w).^3+1.1524e-005*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                         -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(y(:,2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(y(:,2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(y(:,2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                         +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(y(:,2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                         -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(y(:,2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(y(:,2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(y(:,2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(y(:,2)+w0)./v_w).^4 ...
                         +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(y(:,2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(y(:,2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(y(:,2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(y(:,2)+w0)./v_w).^4))));


phaseout.dynamics = dy_dt;

if strcmp(input.auxdata.objfun,'abs')
    phaseout.integrand = abs(P_stl-P_mech);
elseif strcmp(input.auxdata.objfun,'min')
    phaseout.integrand = min(P_stl,P_mech)-min(0,P_stl-P_mech).*(P_stl-P_mech);
end

end
%---------------------------------%
% END: WindTurbineContinuous_GPOPS.m %
%---------------------------------%

%---------------------------------%
% BEGIN: WindTurbineEndpoint_GPOPS.m %
%---------------------------------%
function output = WindTurbineEndpoint_GPOPS(input)
q = input.phase.integral;

% P_inpf = input.phase.finalstate(5);
% P_inp0 = input.phase.initialstate(5);
% beta = input.parameter;

output.objective = q;
end
%---------------------------------%
% END: WindTurbineEndpoint_GPOPS.m %
%---------------------------------%

