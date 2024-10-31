close all;clear all;
format long;
results_structure_all = [];
profile_types = {'Ramp','Gaussian'};
u0_vals = [4,1.25];
for profile_index = 1:length(profile_types)     
        %Wind parameters
        v0 = 11;
        profile_type = char(profile_types(profile_index));
        std = 0.5;
        wind_param.profile_type = profile_type; %'Ramp','Gaussian'
        wind_param.v0 = v0; %initial wind speed
        wind_param.mu = 20; %mean value for the Gaussian profile
        wind_param.std = std;%standard deviation value for the Gaussian profile
        %slope values for the ramp profile
        wind_param.wind_speed_initial = wind_param.v0;
        wind_param.wind_speed_final = wind_param.v0+2;
        wind_param.wind_speed_initial_time = 19;
        wind_param.wind_speed_final_time = 21;
        wind_param.ramp_slope = (wind_param.wind_speed_final-wind_param.wind_speed_initial)/(wind_param.wind_speed_final_time-wind_param.wind_speed_initial_time);
        xo = 1*ones(11,1);
        options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state

        u0 = u0_vals(profile_index);
        g= @(x0) steady_states(x0,u0,wind_param); 
        [x0,~] = fsolve(g,xo,options_steady_states); % finding equilibrium
        [u_opt_all_LD,grad_all_LD,fval_all_LD,topt_LD,uopt_LD,xopt_LD,x0_LD,P_elec_opt_LD,P_mech_opt_LD,v_wind_LD] = SequentialWindTurbineDAE(v0,wind_param,std,u0,x0,'on');

        results_structure.u0_LD = u0;
        results_structure.u_opt_all_LD = u_opt_all_LD;
        results_structure.grad_all_LD = grad_all_LD;
        results_structure.fval_all_LD = fval_all_LD;
        results_structure.topt_LD = topt_LD;
        results_structure.uopt_LD = uopt_LD;
        results_structure.xopt_LD = xopt_LD;
        results_structure.P_mech_opt_LD = P_mech_opt_LD;
        results_structure.P_elec_opt_LD = P_elec_opt_LD;
        
        results_structure.wind_param = wind_param;
        results_structure.x0 = x0;
        results_structure_all = [results_structure_all,results_structure];
    
end

function f = steady_states(x_guess,u0,wind_param)

    u=u0;
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
    %******************************
    R1=0.02;
    X1=0.0243;
    Xtr=0.00557;
    E=1.0164;
    R=R1;
    X=X1+Xtr;

    V_w = wind_param.v0;
    v_w=V_w;

    A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
    B = -(2.*x_guess(10).*R + 2*X.*x_guess(9)./Xeq + 2*(R^2+X^2).*x_guess(9)./Xeq);
    C = (R^2+X^2)/Xeq + (R^2+X^2).*x_guess(10).^2 - E^2;
    V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);
       
    f=[ ...
        ((-x_guess(11)*x_guess(10))/(x_guess(1)+1)-x_guess(3)*Ktg-Dtg*(x_guess(1)-x_guess(2)))/(2*Hg);    
    ((Cmech*V_w^3*((-4.1909e-001+2.1808e-001*(k_b*(x_guess(2)+w0)/v_w)-1.2406e-002*(k_b*(x_guess(2)+w0)/v_w)^2-1.3365e-004*(k_b*(x_guess(2)+w0)/v_w)^3+1.1524e-005*(k_b*(x_guess(2)+w0)/v_w)^4 ...
                                  -6.7606e-002*u    +6.0405e-002*u*(k_b*(x_guess(2)+w0)/v_w)     -1.3934e-002*u*(k_b*(x_guess(2)+w0)/v_w)^2    +1.0683e-003*u*(k_b*(x_guess(2)+w0)/v_w)^3    -2.3895e-005*u*(k_b*(x_guess(2)+w0)/v_w)^4 ...
                                  +1.5727e-002*u^2  -1.0996e-002*u^2*(k_b*(x_guess(2)+w0)/v_w)   +2.1495e-003*u^2*(k_b*(x_guess(2)+w0)/v_w)^2  -1.4855e-004*u^2*(k_b*(x_guess(2)+w0)/v_w)^3  +2.7937e-006*u^2*(k_b*(x_guess(2)+w0)/v_w)^4 ...
                                  -8.6018e-004*u^3  +5.7051e-004*u^3*(k_b*(x_guess(2)+w0)/v_w)   -1.0479e-004*u^3*(k_b*(x_guess(2)+w0)/v_w)^2  +5.9924e-006*u^3*(k_b*(x_guess(2)+w0)/v_w)^3  -8.9194e-008*u^3*(k_b*(x_guess(2)+w0)/v_w)^4 ...
                                  +1.4787e-005*u^4  -9.4839e-006*u^4*(k_b*(x_guess(2)+w0)/v_w)   +1.6167e-006*u^4*(k_b*(x_guess(2)+w0)/v_w)^2  -7.1535e-008*u^4*(k_b*(x_guess(2)+w0)/v_w)^3  +4.9686e-010*u^4*(k_b*(x_guess(2)+w0)/v_w)^4))...
              /(x_guess(2)+w0))+D_tg*(x_guess(1)-x_guess(2))+k_tg*x_guess(3))/(2*H);
    wb*(x_guess(1)-x_guess(2));
    x_guess(1)+1-1.2;
    ((x_guess(1)+1)*(Kptrq*(x_guess(1)+1-1.2) + Kitrq*x_guess(4))-x_guess(5))/Tpc;
    (x_guess(11)*x_guess(8)-x_guess(6))/Tpwr;
    (tan(PFE_ref)*x_guess(6)-(x_guess(11))*(x_guess(9)-x_guess(11))/Xeq)*Kqi;
    (x_guess(7)-x_guess(11))*Kvi;
    (x_guess(8)-x_guess(9))/0.02;
    ((x_guess(5)/x_guess(11))-x_guess(10))/0.02;
    (x_guess(11))^4-(2*(x_guess(11)*x_guess(10)*R+x_guess(11)*(x_guess(9)-x_guess(11))/Xeq*X)+E^2)*(x_guess(11))^2+(R^2+X^2)*((x_guess(11)*x_guess(10))^2+(x_guess(11)*(x_guess(9)-x_guess(11))/Xeq)^2)];

end
function [u_opt_all,grad_all,fval_all,topt,uopt,xopt,x0,P_elec_opt,P_mech_opt,v_wind] = SequentialWindTurbineDAE(v0,wind_param,std,u0,x0,LD_in)

format long;

x0 = [x0(1:10);0;x0(11)]; %add the intial condition for x11, x11(0) = 0
xf=x0;
LD_set = {LD_in};

step_size = 0.2;
Algorithm = 'sqp';
history = [];
gradient_numeric = [];
fval_obj = [];

% Time Horizon and Initial State
t0 = 18; % Initial time
tf = 22; % Final time
xmin = [0.1;0.1;-2;-3;0.7;0.7;0.9;0.5;0.5;0.7];
xmax = [0.4;0.3;2;3;1.5;1.3;1.1;1.45;1.45;1.45];
ns = (tf-t0)/step_size;
M_p = eye(ns); %directions matrix for LD derivatives
ts = t0:(tf-t0)/ns:tf; % Time stages (equipartition)
% Initial Guess and Bounds for the Parameters
u0 = u0*ones(ns,1);
uL = 0*ones(ns,1);
uU = 26*ones(ns,1);
M_DAE=eye(12+12*ns);
for i=[12:12:(ns*12)]
    M_DAE(i,i)=0;
end
% Options for ODE & NLP Solvers
LD = char(LD_set(1));
if strcmp(LD,'off')
    GradObj = 'off';
    linecolor = 'r';
    optODE=odeset('Mass',M_DAE(1:12,1:12),'RelTol',1e-10);
else
    GradObj = 'on';
    linecolor = 'k';
    optODE=odeset('Mass',M_DAE,'RelTol',1e-10);
end

optNLP = optimset( 'Algorithm',Algorithm,'Hessian','bfgs','LargeScale', 'off', 'GradObj', GradObj, 'GradConstr', 'off',...
    'DerivativeCheck', 'off', 'Display', 'iter-detailed', 'TolX', 1e-9,...
    'TolFun', 1e-9, 'TolCon', 1e-6, 'MaxFunEval', 30000, 'Maxiter', 1e+03,'OutputFcn',@myoutput );

% Sequential Approach of Dynamic Optimization
[ uopt ] = fmincon( @(us)obj(x0,ns,ts,us,optODE, wind_param,LD,M_p), u0, [], [], [], [],...
    uL, uU, @(us)ctr(x0,ns,ts,us,optODE,xmin,xmax,xf, wind_param,LD, M_p), optNLP);

plots( x0, ns, ts, uopt, optODE, wind_param, Algorithm, GradObj, LD, linecolor);




    function [ J, dJ ] = obj( x0, ns, ts, us, optODE, wind_param, LD, M_p )
        if strcmp(LD,'off') %solver calculates gradient automatically
            f = fun( x0, ns, ts, us, optODE, wind_param, LD, M_p );
            J = f(11);
        else %provide gradient to the solver
            [f,df] = fun( x0, ns, ts, us, optODE, wind_param, LD, M_p );
            J = f(11);
            dJ = df(11,:)';
        end
    end


    function [ c, ceq, dc, dceq ] = ctr( x0, ns, ts, us, optODE, xmin, xmax, xf, wind_param,LD, M_p) %no constraints
        if nargout == 2
            f = fun( x0, ns, ts, us, optODE, wind_param,LD, M_p );
            ceq = [];
            c = [];
        else
            [f,df] = fun( x0, ns, ts, us, optODE, wind_param,LD, M_p );
            ceq = [];
            dceq = [];
            c = [];
            dc = [];
        end
    end

    function [ f, df ] = fun( x0, ns, ts, us, optODE, wind_param, LD, M_p )
        % Calculate function values only
        if strcmp(LD,'off')
            % Forward state integration
            z0 = x0;
            for ks = 1:ns
                ode = @(t,x)windturbine_dynamics(t,x,us,ks,ts,wind_param);
                [tspan,zs] = ode15s( ode, [ts(ks),ts(ks+1)], z0, optODE );
                z0 = zs(end,:)';
            end
            f = zs(end,:)'; %get the objective function J = x11(tf)
            
            % Calculate both function and gradient values
        else
            % Forward state & sensitivity integration
            z0 = [ x0; zeros(12*ns,1) ];
            for ks = 1:ns
                ode = @(t,x)windturbine_dynamics_and_sensentivities(t,x,us,ks,ts,wind_param,LD, M_p);
                [tspan,zs] = ode15s( ode, [ts(ks),ts(ks+1)], z0, optODE );
                z0 = zs(end,:)';
            end
            % Functions & Gradients
            f = zs(end,1:12)'; %get the objective function J = x11(tf)
            df = [];
            for is = 1:ns %collect senstivities
                df(1:12,is) = zs(end,is*12+1:is*12+12);
            end
        end
    end


    function plots( x0, ns, ts, us, optODE,wind_param, Algorithm, GradObj, LD, linecolor)
        % Forward state integration; store optimal state & control
        z0 = x0;
        topt = [];
        xopt = [];
        uopt = [];
        P_elec_opt = [];
        P_mech_opt = [];
        
        
        v_wind = wind_param.v0;
        Xeq=0.8;
        R1=0.02;
        X1=0.0243;
        Xtr=0.00557;
        E=1.0164;
        R=R1;
        X=X1+Xtr;
        
        v_step = wind_param.wind_speed_initial;
        optODE=odeset('Mass',M_DAE(1:12,1:12),'RelTol',1e-10);
        for ks = 1:ns
            ode = @(t,x)windturbine_dynamics(t,x,us,ks,ts,wind_param);
            [tspan,zs] = ode15s( ode, [ts(ks),ts(ks+1)], z0, optODE );
            z0 = zs(end,:)';
            for tt = 1:size(tspan,1)
                topt = [ topt; tspan(tt) ];
                xopt = [ xopt; zs(tt,:) ];
                uopt = [ uopt; us(ks)];
                if strcmp(wind_param.profile_type,'Gaussian')
                    rand_profile = (1/(wind_param.std*sqrt(2*pi)))*exp(-1/2*((tspan(tt)-wind_param.mu)/wind_param.std).^2);
                    v_step = wind_param.wind_speed_initial + rand_profile;
                    v_wind = [v_wind; wind_param.wind_speed_initial + rand_profile];
                else
                    if (tspan(tt) < wind_param.wind_speed_initial_time)
                        ramp_profile = 0;
                    elseif (tspan(tt) > wind_param.wind_speed_final_time)
                        ramp_profile = wind_param.ramp_slope*(wind_param.wind_speed_final_time-wind_param.wind_speed_initial_time);
                    else
                        ramp_profile = wind_param.ramp_slope*(tspan(tt)-wind_param.wind_speed_initial_time);
                    end
                    v_step = wind_param.wind_speed_initial + ramp_profile;
                    v_wind = [v_wind; v_step];
                end
                
                P_elec_step = zs(tt,6);
                
                P_elec_opt = [ P_elec_opt; P_elec_step ];
                
                v_w = v_step;
                
                w0=1;
                D_tg=1.5;
                k_tg=1.11;
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
                R1=0.02;
                X1=0.0243;
                Xtr=0.00557;
                E=1.0164;
                R=R1;
                X=X1+Xtr;
                
                P_mech_step = Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(zs(tt,2)+w0)./v_w)-1.2406e-002*(k_b.*(zs(tt,2)+w0)./v_w).^2-1.3365e-004*(k_b.*(zs(tt,2)+w0)./v_w).^3+1.1524e-005*(k_b.*(zs(tt,2)+w0)./v_w).^4 ...
                    -6.7606e-002.*us(ks) +6.0405e-002.*us(ks).*(k_b.*(zs(tt,2)+w0)./v_w)     -1.3934e-002.*us(ks).*(k_b.*(zs(tt,2)+w0)./v_w).^2    +1.0683e-003.*us(ks).*(k_b.*(zs(tt,2)+w0)./v_w).^3    -2.3895e-005.*us(ks).*(k_b.*(zs(tt,2)+w0)./v_w).^4 ...
                    +1.5727e-002.*us(ks).^2  -1.0996e-002.*us(ks).^2.*(k_b.*(zs(tt,2)+w0)./v_w)   +2.1495e-003.*us(ks).^2.*(k_b.*(zs(tt,2)+w0)./v_w).^2  -1.4855e-004.*us(ks).^2.*(k_b.*(zs(tt,2)+w0)./v_w).^3  +2.7937e-006.*us(ks).^2.*(k_b.*(zs(tt,2)+w0)./v_w).^4 ...
                    -8.6018e-004.*us(ks).^3  +5.7051e-004.*us(ks).^3.*(k_b.*(zs(tt,2)+w0)./v_w)   -1.0479e-004.*us(ks).^3.*(k_b.*(zs(tt,2)+w0)./v_w).^2  +5.9924e-006.*us(ks).^3.*(k_b*(zs(tt,2)+w0)./v_w).^3  -8.9194e-008.*us(ks).^3.*(k_b.*(zs(tt,2)+w0)./v_w).^4 ...
                    +1.4787e-005.*us(ks).^4  -9.4839e-006.*us(ks).^4.*(k_b.*(zs(tt,2)+w0)./v_w)   +1.6167e-006.*us(ks).^4.*(k_b.*(zs(tt,2)+w0)./v_w).^2  -7.1535e-008.*us(ks).^4.*(k_b*(zs(tt,2)+w0)./v_w).^3  +4.9686e-010.*us(ks).^4.*(k_b.*(zs(tt,2)+w0)./v_w).^4));
                
                P_mech_opt = [P_mech_opt;P_mech_step];
            end
        end
        v_wind(112) = [];
        
        font_size = 20;
        figure('Name','Pmech','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
        hold on
        plot(topt,P_mech_opt,linecolor)
        hold on
        % plot(topt,P_elec_opt,'k--')
        xlabel('t')
        ylabel('Pmech')
        
        
        figure('Name','u(t)','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
        hold on
        plot(topt,uopt,linecolor)
        hold on
        xlabel('t')
        ylabel('u(t)')
        
        
        % smoothing_param = 0.9;
        % uopt_smooth = csaps(topt,uopt,smoothing_param);
        % figure(u_smooth_fig)
        % hold on
        % fnplt(uopt_smooth,linecolor)
        % hold on
        % xlabel('t')
        % ylabel('u(t) interpolated')
        
        
        figure('Name','v','DefaultAxesFontSize',font_size,'defaultLineLineWidth',3,'DefaultAxesTitleFontWeight','bold');
        plot(topt,v_wind,linecolor)
        xlabel('t')
        ylabel('v_w')
        
    end

    function [dx_dt] = windturbine_dynamics(t,x,us,ks,ts,wind_param)
        
        %Get the current input u
        u = us(ks);
        ns = length(us);
        %for wind_param.mu extra state equal to the integrand in the objective functional to transform from Lagrange to Mayer
        Xeq=0.8;
        R1=0.02;
        X1=0.0243;
        Xtr=0.00557;
        E=1.0164;
        R=R1;
        X=X1+Xtr;
        A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
        B = -(2.*x(10).*R + 2*X.*x(9)./Xeq + 2*(R^2+X^2).*x(9)./Xeq);
        C = (R^2+X^2)/Xeq + (R^2+X^2).*x(10,:).^2 - E^2;
        V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);
        P_inp = x(5);
        
        if strcmp(wind_param.profile_type,'Gaussian')
            rand_profile = (1/(wind_param.std*sqrt(2*pi)))*exp(-1/2*((t-wind_param.mu)/wind_param.std).^2);
            wind_speed_profile = wind_param.v0 + rand_profile;
        else
            wind_speed_profile = wind_param.wind_speed_initial;
            if t < wind_param.wind_speed_initial_time
                ramp_profile = 0;
            elseif t > wind_param.wind_speed_final_time
                ramp_profile = wind_param.ramp_slope*(wind_param.wind_speed_final_time -wind_param.wind_speed_initial_time);
            else
                ramp_profile = wind_param.ramp_slope*(t-wind_param.wind_speed_initial_time);
            end
            wind_speed_profile = wind_speed_profile + ramp_profile;
        end
        v_w = wind_speed_profile';
        
        w0=1;
        D_tg=1.5;
        k_tg=1.11;
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
        R1=0.02;
        X1=0.0243;
        Xtr=0.00557;
        E=1.0164;
        R=R1;
        X=X1+Xtr;
        Pelec = x(12)*x(10);
        if v_w <6
            w_ref = -0.75*Pelec^2+1.59*Pelec+0.63;
        else
            w_ref = 1.2;
        end
        P_stl = 1;
        P_mech = Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(x(2)+w0)./v_w)-1.2406e-002*(k_b.*(x(2)+w0)./v_w).^2-1.3365e-004*(k_b.*(x(2)+w0)./v_w).^3+1.1524e-005*(k_b.*(x(2)+w0)./v_w).^4 ...
            -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(x(2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(x(2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(x(2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(x(2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(x(2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(x(2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(x(2)+w0)./v_w).^4 ...
            -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(x(2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(x(2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(x(2)+w0)./v_w).^4));
        
        %Begin of dynamics
        dx_dt(1)=(-(x(12)*x(10)./(x(1)+1))-x(3).*Ktg-Dtg.*(x(1)-x(2)))./(2*Hg);
        dx_dt(2)=((Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(x(2)+w0)./v_w)-1.2406e-002*(k_b.*(x(2)+w0)./v_w).^2-1.3365e-004*(k_b.*(x(2)+w0)./v_w).^3+1.1524e-005*(k_b.*(x(2)+w0)./v_w).^4 ...
            -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(x(2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(x(2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(x(2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(x(2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(x(2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(x(2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(x(2)+w0)./v_w).^4 ...
            -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(x(2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(x(2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(x(2)+w0)./v_w).^4))...
            ./(x(2)+w0))+D_tg.*(x(1)-x(2))+k_tg*x(3))./(2*H);
        dx_dt(3)=wb.*(x(1)-x(2));
        dx_dt(4)=x(1)+1-w_ref;
        dx_dt(5)=((x(1)+1).*(Kptrq.*(x(1)+1-w_ref) + Kitrq.*x(4))-x(5))./Tpc;
        dx_dt(6)=(x(12)*x(10)-x(6))/Tpwr;
        dx_dt(7)=(tan(PFE_ref)*x(6)-(x(12)*(x(9)-x(12))/Xeq))*Kqi;
        dx_dt(8)=(x(7)-x(12))*Kvi;
        dx_dt(9)=(x(8)-x(9))./0.02;
        dx_dt(10)=((x(5)./x(12))-x(10))/0.02;
        x11_dot = -P_stl+(P_mech - P_stl)^2;
        dx_dt(11) = x11_dot;
        dx_dt(12)=(x(12))^4-(2*(x(12)*x(10)*R+x(12)*(x(9)-x(12))/Xeq*X)+E^2)*(x(12))^2+(R^2+X^2)*((x(12)*x(10))^2+(x(12)*(x(9)-x(12))/Xeq)^2);
        %End of dynamics
        
        dx_dt = dx_dt';
        
    end

    function [dx_dt] = windturbine_dynamics_and_sensentivities(t,x,us,ks,ts,wind_param, LD, M_p)
        
        %Get the current input u
        u = us(ks);
        ns = length(us);
        %for wind_param.mu extra state equal to the integrand in the objective functional to transform from Lagrange to Mayer
        Xeq=0.8;
        R1=0.02;
        X1=0.0243;
        Xtr=0.00557;
        E=1.0164;
        R=R1;
        X=X1+Xtr;
        A = 1 + 2*X/Xeq + (R^2+X^2)/Xeq;
        B = -(2.*x(10).*R + 2*X.*x(9)./Xeq + 2*(R^2+X^2).*x(9)./Xeq);
        C = (R^2+X^2)/Xeq + (R^2+X^2).*x(10,:).^2 - E^2;
        V = (-1.*B+sqrt(B.^2-4*A*C))./(2*A);
        
        P_inp = x(5); %same as Pelec
        
        if strcmp(wind_param.profile_type,'Gaussian')
            rand_profile = (1/(wind_param.std*sqrt(2*pi)))*exp(-1/2*((t-wind_param.mu)/wind_param.std).^2);
            wind_speed_profile = wind_param.v0 + rand_profile;
        else %Ramp wind profile
            wind_speed_profile = wind_param.wind_speed_initial;
            if t < wind_param.wind_speed_initial_time
                ramp_profile = 0;
            elseif t > wind_param.wind_speed_final_time
                ramp_profile = wind_param.ramp_slope*(wind_param.wind_speed_final_time -wind_param.wind_speed_initial_time);
            else
                ramp_profile = wind_param.ramp_slope*(t-wind_param.wind_speed_initial_time);
            end
            wind_speed_profile = wind_speed_profile + ramp_profile;
        end
        v_w = wind_speed_profile';
        
        w0=1;
        D_tg=1.5;
        k_tg=1.11;
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
        R1=0.02;
        X1=0.0243;
        Xtr=0.00557;
        E=1.0164;
        R=R1;
        X=X1+Xtr;
        
        P_stl = 1;
        P_mech = Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(x(2)+w0)./v_w)-1.2406e-002*(k_b.*(x(2)+w0)./v_w).^2-1.3365e-004*(k_b.*(x(2)+w0)./v_w).^3+1.1524e-005*(k_b.*(x(2)+w0)./v_w).^4 ...
            -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(x(2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(x(2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(x(2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(x(2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(x(2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(x(2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(x(2)+w0)./v_w).^4 ...
            -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(x(2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(x(2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(x(2)+w0)./v_w).^4));
        
        %Begin of dynamics
        dx_dt(1)=(-(x(12)*x(10)./(x(1)+1))-x(3).*Ktg-Dtg.*(x(1)-x(2)))./(2*Hg);
        dx_dt(2)=((Cmech*v_w.^3.*((-4.1909e-001+2.1808e-001*(k_b.*(x(2)+w0)./v_w)-1.2406e-002*(k_b.*(x(2)+w0)./v_w).^2-1.3365e-004*(k_b.*(x(2)+w0)./v_w).^3+1.1524e-005*(k_b.*(x(2)+w0)./v_w).^4 ...
            -6.7606e-002.*u +6.0405e-002.*u.*(k_b.*(x(2)+w0)./v_w)     -1.3934e-002.*u.*(k_b.*(x(2)+w0)./v_w).^2    +1.0683e-003.*u.*(k_b.*(x(2)+w0)./v_w).^3    -2.3895e-005.*u.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.5727e-002.*u.^2  -1.0996e-002.*u.^2.*(k_b.*(x(2)+w0)./v_w)   +2.1495e-003.*u.^2.*(k_b.*(x(2)+w0)./v_w).^2  -1.4855e-004.*u.^2.*(k_b.*(x(2)+w0)./v_w).^3  +2.7937e-006.*u.^2.*(k_b.*(x(2)+w0)./v_w).^4 ...
            -8.6018e-004.*u.^3  +5.7051e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w)   -1.0479e-004.*u.^3.*(k_b.*(x(2)+w0)./v_w).^2  +5.9924e-006.*u.^3.*(k_b*(x(2)+w0)./v_w).^3  -8.9194e-008.*u.^3.*(k_b.*(x(2)+w0)./v_w).^4 ...
            +1.4787e-005.*u.^4  -9.4839e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w)   +1.6167e-006.*u.^4.*(k_b.*(x(2)+w0)./v_w).^2  -7.1535e-008.*u.^4.*(k_b*(x(2)+w0)./v_w).^3  +4.9686e-010.*u.^4.*(k_b.*(x(2)+w0)./v_w).^4))...
            ./(x(2)+w0))+D_tg.*(x(1)-x(2))+k_tg*x(3))./(2*H);
        dx_dt(3)=wb.*(x(1)-x(2));
        dx_dt(4)=x(1)+1-1.2;
        dx_dt(5)=((x(1)+1).*(Kptrq.*(x(1)+1-1.2) + Kitrq.*x(4))-x(5))./Tpc;
        dx_dt(6)=(x(12)*x(10)-x(6))/Tpwr;
        dx_dt(7)=(tan(PFE_ref)*x(6)-(x(12)*(x(9)-x(12))/Xeq))*Kqi;
        dx_dt(8)=(x(7)-x(12))*Kvi;
        dx_dt(9)=(x(8)-x(9))./0.02;
        dx_dt(10)=((x(5)./x(12))-x(10))/0.02;
        x11_dot = -(min(P_stl,P_mech)) + min(0,P_stl-P_mech)*(P_stl-P_mech);
        dx_dt(11) = x11_dot;
        dx_dt(12)=(x(12))^4-(2*(x(12)*x(10)*R+x(12)*(x(9)-x(12))/Xeq*X)+E^2)*(x(12))^2+(R^2+X^2)*((x(12)*x(10))^2+(x(12)*(x(9)-x(12))/Xeq)^2);
        %End of dynamics
        
        dx_dt = dx_dt';
        
        % Append sensitivity system
        for is = 1:ns
            if is == ks
                uws = 1;
            else
                uws = 0;
            end
            x1 = x(1);x2 = x(2);x3 = x(3);x4 = x(4);x5 = x(5);x6 = x(6);x7 = x(7);x8 = x(8);x9 = x(9);x10 = x(10);x11 = x(11);x12 = x(12);
            x1u = x(12*is+1);x2u = x(12*is+2);x3u = x(12*is+3);x4u = x(12*is+4);x5u = x(12*is+5);x6u = x(12*is+6);x7u = x(12*is+7);x8u = x(12*is+8);
            x9u = x(12*is+9);x10u = x(12*is+10);x11u = x(12*is+11);x12u = x(12*is+12);
            Sx = [x1u;x2u;x3u;x4u;x5u;x6u;x7u;x8u;x9u;x10u;x11u;x12u];
            
            dPmech_dx = [ 0, Cmech*v_w^3*((0.2181*k_b)/v_w - (0.0124*k_b^2*(2*w0 + 2*x2))/v_w^2 + (0.0604*k_b*u)/v_w - (0.0110*k_b*u^2)/v_w + (5.7051e-04*k_b*u^3)/v_w - (9.4839e-06*k_b*u^4)/v_w - (4.0095e-04*k_b^3*(w0 + x2)^2)/v_w^3 + (4.6096e-05*k_b^4*(w0 + x2)^3)/v_w^4 + (0.0032*k_b^3*u*(w0 + x2)^2)/v_w^3 - (9.5580e-05*k_b^4*u*(w0 + x2)^3)/v_w^4 - (0.0139*k_b^2*u*(2*w0 + 2*x2))/v_w^2 - (4.4565e-04*k_b^3*u^2*(w0 + x2)^2)/v_w^3 + (1.7977e-05*k_b^3*u^3*(w0 + x2)^2)/v_w^3 - (2.1460e-07*k_b^3*u^4*(w0 + x2)^2)/v_w^3 + (1.1175e-05*k_b^4*u^2*(w0 + x2)^3)/v_w^4 - (3.5678e-07*k_b^4*u^3*(w0 + x2)^3)/v_w^4 + (1.9874e-09*k_b^4*u^4*(w0 + x2)^3)/v_w^4 + (0.0021*k_b^2*u^2*(2*w0 + 2*x2))/v_w^2 - (1.0479e-04*k_b^2*u^3*(2*w0 + 2*x2))/v_w^2 + (1.6167e-06*k_b^2*u^4*(2*w0 + 2*x2))/v_w^2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            dPmech_du = Cmech*v_w^3*(0.0315*u - 0.0026*u^2 + 5.9148e-05*u^3 - (0.0139*k_b^2*(w0 + x2)^2)/v_w^2 + (0.0011*k_b^3*(w0 + x2)^3)/v_w^3 - (2.3895e-05*k_b^4*(w0 + x2)^4)/v_w^4 + (0.0604*k_b*(w0 + x2))/v_w + (0.0043*k_b^2*u*(w0 + x2)^2)/v_w^2 - (2.9710e-04*k_b^3*u*(w0 + x2)^3)/v_w^3 + (5.5874e-06*k_b^4*u*(w0 + x2)^4)/v_w^4 - (0.0220*k_b*u*(w0 + x2))/v_w - (3.1437e-04*k_b^2*u^2*(w0 + x2)^2)/v_w^2 + (6.4668e-06*k_b^2*u^3*(w0 + x2)^2)/v_w^2 + (1.7977e-05*k_b^3*u^2*(w0 + x2)^3)/v_w^3 - (2.8614e-07*k_b^3*u^3*(w0 + x2)^3)/v_w^3 - (2.6758e-07*k_b^4*u^2*(w0 + x2)^4)/v_w^4 + (1.9874e-09*k_b^4*u^3*(w0 + x2)^4)/v_w^4 + (0.0017*k_b*u^2*(w0 + x2))/v_w - (3.7936e-05*k_b*u^3*(w0 + x2))/v_w - 0.0676);
            
            
            
            dfdx = [             -(0.5000*(Dtg - (x10*x12)/(x1 + 1)^2))/Hg,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (0.5000*Dtg)/Hg, -(0.5000*Ktg)/Hg,                    0,      0,                0,   0,  0,              0, -(0.5000*x12)/(Hg*(x1 + 1)), 0,    -(0.5000*x10)/(Hg*(x1 + 1));
                (0.5000*D_tg)/H, (0.5000*((Cmech*v_w^3*((0.2181*k_b)/v_w - (0.0124*k_b^2*(2*w0 + 2*x2))/v_w^2 + (0.0604*k_b*u)/v_w - (0.0110*k_b*u^2)/v_w + (5.7051e-04*k_b*u^3)/v_w - (9.4839e-06*k_b*u^4)/v_w - (4.0095e-04*k_b^3*(w0 + x2)^2)/v_w^3 + (4.6096e-05*k_b^4*(w0 + x2)^3)/v_w^4 + (0.0032*k_b^3*u*(w0 + x2)^2)/v_w^3 - (9.5580e-05*k_b^4*u*(w0 + x2)^3)/v_w^4 - (0.0139*k_b^2*u*(2*w0 + 2*x2))/v_w^2 - (4.4565e-04*k_b^3*u^2*(w0 + x2)^2)/v_w^3 + (1.7977e-05*k_b^3*u^3*(w0 + x2)^2)/v_w^3 - (2.1460e-07*k_b^3*u^4*(w0 + x2)^2)/v_w^3 + (1.1175e-05*k_b^4*u^2*(w0 + x2)^3)/v_w^4 - (3.5678e-07*k_b^4*u^3*(w0 + x2)^3)/v_w^4 + (1.9874e-09*k_b^4*u^4*(w0 + x2)^3)/v_w^4 + (0.0021*k_b^2*u^2*(2*w0 + 2*x2))/v_w^2 - (1.0479e-04*k_b^2*u^3*(2*w0 + 2*x2))/v_w^2 + (1.6167e-06*k_b^2*u^4*(2*w0 + 2*x2))/v_w^2))/(w0 + x2) - D_tg + (Cmech*v_w^3*(0.0676*u - 0.0157*u^2 + 8.6018e-04*u^3 - 1.4787e-05*u^4 + (0.0124*k_b^2*(w0 + x2)^2)/v_w^2 + (1.3365e-04*k_b^3*(w0 + x2)^3)/v_w^3 - (1.1524e-05*k_b^4*(w0 + x2)^4)/v_w^4 - (0.2181*k_b*(w0 + x2))/v_w + (0.0139*k_b^2*u*(w0 + x2)^2)/v_w^2 - (0.0011*k_b^3*u*(w0 + x2)^3)/v_w^3 + (2.3895e-05*k_b^4*u*(w0 + x2)^4)/v_w^4 - (0.0604*k_b*u*(w0 + x2))/v_w - (0.0021*k_b^2*u^2*(w0 + x2)^2)/v_w^2 + (1.0479e-04*k_b^2*u^3*(w0 + x2)^2)/v_w^2 - (1.6167e-06*k_b^2*u^4*(w0 + x2)^2)/v_w^2 + (1.4855e-04*k_b^3*u^2*(w0 + x2)^3)/v_w^3 - (5.9924e-06*k_b^3*u^3*(w0 + x2)^3)/v_w^3 + (7.1535e-08*k_b^3*u^4*(w0 + x2)^3)/v_w^3 - (2.7937e-06*k_b^4*u^2*(w0 + x2)^4)/v_w^4 + (8.9194e-08*k_b^4*u^3*(w0 + x2)^4)/v_w^4 - (4.9686e-10*k_b^4*u^4*(w0 + x2)^4)/v_w^4 + (0.0110*k_b*u^2*(w0 + x2))/v_w - (5.7051e-04*k_b*u^3*(w0 + x2))/v_w + (9.4839e-06*k_b*u^4*(w0 + x2))/v_w + 0.4191))/(w0 + x2)^2))/H,  (0.5000*k_tg)/H,                    0,      0,                0,   0,  0,              0,                           0, 0,                              0;
                wb,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    -wb,                0,                    0,      0,                0,   0,  0,              0,                           0, 0,                              0;
                1,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0,      0,                0,   0,  0,              0,                           0, 0,                              0;
                (Kitrq*x4 + Kptrq*(x1 + 1) + Kptrq*(x1 - 0.2000))/Tpc,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0, (Kitrq*(x1 + 1))/Tpc, -1/Tpc,                0,   0,  0,              0,                           0, 0,                              0;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0,      0,          -1/Tpwr,   0,  0,              0,                    x12/Tpwr, 0,                       x10/Tpwr;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0,      0, Kqi*tan(PFE_ref),   0,  0, -(Kqi*x12)/Xeq,                           0, 0, Kqi*(x12/Xeq - (x9 - x12)/Xeq);
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0,      0,                0, Kvi,  0,              0,                           0, 0,                           -Kvi;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0,      0,                0,   0, 50,            -50,                           0, 0,                              0;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      0,                0,                    0, 50/x12,                0,   0,  0,              0,                         -50, 0,                 -(50*x5)/x12^2];
            
            dfdu = [ 0;...
                (0.5000*Cmech*v_w^3*(0.0315*u - 0.0026*u^2 + 5.9148e-05*u^3 - (0.0139*k_b^2*(w0 + x2)^2)/v_w^2 + (0.0011*k_b^3*(w0 + x2)^3)/v_w^3 - (2.3895e-05*k_b^4*(w0 + x2)^4)/v_w^4 + (0.0604*k_b*(w0 + x2))/v_w + (0.0043*k_b^2*u*(w0 + x2)^2)/v_w^2 - (2.9710e-04*k_b^3*u*(w0 + x2)^3)/v_w^3 + (5.5874e-06*k_b^4*u*(w0 + x2)^4)/v_w^4 - (0.0220*k_b*u*(w0 + x2))/v_w - (3.1437e-04*k_b^2*u^2*(w0 + x2)^2)/v_w^2 + (6.4668e-06*k_b^2*u^3*(w0 + x2)^2)/v_w^2 + (1.7977e-05*k_b^3*u^2*(w0 + x2)^3)/v_w^3 - (2.8614e-07*k_b^3*u^3*(w0 + x2)^3)/v_w^3 - (2.6758e-07*k_b^4*u^2*(w0 + x2)^4)/v_w^4 + (1.9874e-09*k_b^4*u^3*(w0 + x2)^4)/v_w^4 + (0.0017*k_b*u^2*(w0 + x2))/v_w - (3.7936e-05*k_b*u^3*(w0 + x2))/v_w - 0.0676))/(H*(w0 + x2));...
                0;...
                0;...
                0;...
                0;...
                0;...
                0;...
                0;...
                0 ...
                ];
            Sx_dot=zeros(11,1);
            Sx_dot(1:10) = dfdu*uws + dfdx*Sx;
            
            %LD derivative for nonsmooth x11, M and N are directions matrices used in the
            %chain rule of LD derivatives
            M = dPmech_du*uws*M_p(ks,is)'; %here M is a scalar R1 * R1 -> R1
            N = dPmech_dx*Sx; %here N is a scalar R1x11 * R11x1 -> R1
            dSx11_dt_LD = -SLmin([P_stl; 0*M + 0*N],[P_mech; M+N])...
                + SLmin([0; 0*M + 0*N],[P_stl-P_mech; -M-N])*(P_stl-P_mech) ...
                - min(0,P_stl-P_mech)*(M+N);
            Sx_dot(11) = dSx11_dt_LD;
            
            dgdu = 0;
            dgdx = [ 0, 0, 0, 0, 0, 0, 0, 0, (x12^2*(2*x9 - 2*x12)*(R^2 + X^2))/Xeq^2 - (2*X*x12^3)/Xeq, 2*x10*x12^2*(R^2 + X^2) - 2*R*x12^3, 0, (R^2 + X^2)*(2*x10^2*x12 + (2*x12*(x9 - x12)^2)/Xeq^2 - (x12^2*(2*x9 - 2*x12))/Xeq^2) - x12^2*(2*R*x10 + (2*X*(x9 - x12))/Xeq - (2*X*x12)/Xeq) - 2*x12*(E^2 + 2*R*x10*x12 + (2*X*x12*(x9 - x12))/Xeq) + 4*x12^3];
            Sy_dot = dgdu*uws + dgdx*Sx;
            dx_dt = [ ...
                dx_dt;
                Sx_dot;
                Sy_dot;
                ];
        end
    end


    function SLmin_output = SLmin(xM1,yM2)
        % Shifted Lmin
        % Returns the minimum of 2 vectors x&y (where x is the first point augmented with the first row of the directions matrix M and y is the second point
        % augmented with the second row of the direction matrix M ,after removing the first element of the vector
        % Vectors x&y must have the same size
        
        SLmin_output = xM1(2:end);
        for k=1:size(xM1,2)
            if xM1(k)<yM2(k)
                SLmin_output = xM1(2:end);
                break
            elseif xM1(k)>yM2(k)
                SLmin_output = yM2(2:end);
                break
            end
        end
    end

    function stop = myoutput(x,optimvalues,state) %This function to get the log of values of interest within each iteration of the NLP solver
        stop = false;
        if state == 'iter'
            history = [history; x];
            gradient_numeric = [gradient_numeric;optimvalues.gradient];
            fval_obj = [fval_obj;optimvalues.fval];
        end
    end
u_opt_all = history;
grad_all = gradient_numeric;
fval_all = fval_obj;
end