function FigS10E
% function my_EN

    clear 
    close all
    clc
    warning('off','all');  

global param

%% Simulation Time
tspan = [0 630];
param.dt = 0.001;     % how many second per data point for the solver and recording 
t = tspan(1):param.dt:tspan(2);

param.T_scale = 0.0998;

%% Simulation Domain
param.Np = 301; % number of spatial point along the perimeter
% param.dx = 1;
param.L  = 30;
x = linspace(0,param.L,param.Np);
param.dx = mean(diff(x));
param.x = x;


%% Parameters
% parameters in terms of the reaction rates

param.a1 = 0.083;
param.a2 = 8.33;
param.a3 = 93.7;
param.a4 = 2.4;
param.a5 = 4; 
param.ub = 0.735;
param.c1 = 0.1;
param.c2 = 2.1750;


param.D_F = 9;
param.D_R = 18;
param.D_B = 4.5;

param.A1 = 1.2;
param.A2 = 0.8;

param.b1 = 1;
param.b2 = 1;
param.b3 = 1;

param.C1 = 4.5;
param.C2 = 8.5;

param.D_W = 0.1;
param.D_Z = 100;



param.phi = 1;


param.D_Cnmf = 10;
param.D_Cnmb = 0.001;

%% Initial Condition

Sol1 = findeqm;
% return
X0 = [Sol1.x;Sol1.y;1/Sol1.x;...
      0;0;...
      0;0];


%% SDE toolbox
rng(100)
timeRange = t;
initvalues = X0;

problem = 'FigS10E'; % name of the sde file
numsim  = param.Np; % number of spatial points
sdetype = 'Ito';
numdepvars = 3+2+2;

% Noise Characteristics
input.baseMean   = 0;    
input.baseSigma  = 0.3;  

tic
Output = SDE_euler(initvalues,problem,timeRange,numdepvars,numsim,sdetype,param,input);
%%
T_sampling = 1/param.dt/10;
F_SDEtoolbox    = Output(1:T_sampling:end,1:numdepvars:end);
R_SDEtoolbox    = Output(1:T_sampling:end,2:numdepvars:end);
B_SDEtoolbox    = Output(1:T_sampling:end,3:numdepvars:end);


W_SDEtoolbox    = Output(1:T_sampling:end,4:numdepvars:end);
Z_SDEtoolbox    = Output(1:T_sampling:end,5:numdepvars:end);

Cnmf_SDEtoolbox = Output(1:T_sampling:end,6:numdepvars:end);
Cnmb_SDEtoolbox = Output(1:T_sampling:end,7:numdepvars:end);

toc

T = t(1:T_sampling:end);

%% Plot

d = 30;
close all
T_initial = 30;



figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
hold on
surf(T-T_initial,x,circshift(B_SDEtoolbox',d,1),'edgecolor','interp')
colorbar
caxis([floor(min(min(B_SDEtoolbox(501:end,:)))) floor(max(max(B_SDEtoolbox(501:end,:))))])
colormap(gray)
xlim([0 tspan(end)-T_initial])
xlabel('time','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title('Figure S10E: B Kymograph','fontsize',20,'fontweight','n')









%%
Cnm_SDEtoolbox = Cnmf_SDEtoolbox + Cnmb_SDEtoolbox;
figure('color','white')
set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
hold on
surf(T-T_initial,x,circshift(Cnm_SDEtoolbox',d,1),'edgecolor','interp')
colorbar
colormap(gray)
xlim([0 tspan(end)-T_initial])
caxis([20 100])
xlabel('Time (A.U.)','fontsize',20)
ylabel('cell perimeter','fontsize',20)
title('Figure S10E: $$C^-_m =  C^-+_{mf}+ C^-_{mb} \quad \mathrm{Kymograph}$$ ','Interpreter','Latex','fontweight','n','fontsize',20)
%%



end
%% Equilibrium and Nullcline Computation


function Sol1 = findeqm

    global param

    syms F R
    f1 = -(param.a1+param.a2*R).*F + ((param.a3*F.^2)./(param.a4^2+F.^2)+param.ub).*(param.a5-F);

    f5 = -(param.a1+param.a2*R).*F + ((1.2*param.a3*F.^2)./(0.8*param.a4^2+F.^2)+param.ub).*(param.a5-F);

    f7 = -4.5*param.c1*R + param.c2*F*(8.5-R);
    f8 = -1.2*4.5*param.c1*R + param.c2*F*(9.1-R);
    
    F_val = linspace(0.1,6,10001);

    
    F1 = (-param.a1.*F_val + ((param.a3*F_val.^2)./(param.a4^2+F_val.^2)+param.ub).*(param.a5-F_val))./(param.a2*F_val);
    F5 = (-param.a1.*F_val + ((param.A1*param.a3*F_val.^2)./(param.A2*param.a4^2+F_val.^2)+param.ub).*(param.a5-F_val))./(param.a2*F_val);
    F7 = param.c2*F_val*8.5./(4.5*param.c1+param.c2*F_val);


    

    TF1 = islocalmin(abs(F1-F7),'MaxNumExtrema',1);
    index1 =  find(TF1);
    Sol1.x = F_val(index1);
    Sol1.y = F7(index1);

    TF2 = islocalmin(abs(F5-F7),'MaxNumExtrema',3);
    index2 =  find(TF2);
    Sol2.x = F_val(index2);
    Sol2.y = F7(index2);
    


end

