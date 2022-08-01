function Fig7FG

    clear 
    close all
    clc
    warning('off','all');  

    global param

%% Simulation Time


    tspan = [0 600];
    param.dt = 0.001;     % how many second per data point for the solver and recording 
    t = tspan(1):param.dt:tspan(2);

    param.T_scale = 0.1109;
    param.T1 = 0.5;
    param.T2 = 0.5;
    param.Inh = 0.0055;
    
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

param.c2 =2.5;


param.D_F = 11.25;
param.D_R = 22.5;
param.D_B = 4.95;

param.A1 = 0.985;
param.A2 = 0.7344;

param.b1 = 1;
param.b2 = 1;
param.b3 = 1;

param.C1 = 4.5;
param.C2 = 8.5;

param.D_W = 0.1;
param.D_Z = 100;

param.phi = 9;


param.D_Cpmf = 10;
param.D_Cpmb = 0.001;

%% Initial Condition

    Sol1 = findeqm;

    X0 = [Sol1.x;Sol1.y;1/Sol1.x;...
          0;0;...
          0;0];


%% SDE toolbox

    rng(100) % fixing the rng for data reproducibility
    timeRange = t;
    initvalues = X0;

    problem = 'Fig7FG'; % name of the sde file
    numsim  = param.Np; % number of spatial points
    sdetype = 'Ito';
    numdepvars = 3+2+2;

% Noise Characteristics
    input.baseMean   = 0;    
    input.baseSigma  = 0.3;  

    tic
    Output = SDE_euler(initvalues,problem,timeRange,numdepvars,numsim,sdetype,param,input);


    T_sampling = 1/param.dt/10;

    F_SDEtoolbox    = Output(1:T_sampling:end,1:numdepvars:end);
    R_SDEtoolbox    = Output(1:T_sampling:end,2:numdepvars:end);
    B_SDEtoolbox    = Output(1:T_sampling:end,3:numdepvars:end);

    W_SDEtoolbox    = Output(1:T_sampling:end,4:numdepvars:end);
    Z_SDEtoolbox    = Output(1:T_sampling:end,5:numdepvars:end);

    Cpmf_SDEtoolbox = Output(1:T_sampling:end,6:numdepvars:end);
    Cpmb_SDEtoolbox = Output(1:T_sampling:end,7:numdepvars:end);

    toc

    T = t(1:T_sampling:end);

%% Plot


    d = 0;

    T_initial = 0;

% F kymograph    
    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    surf(T-T_initial,x,circshift(F_SDEtoolbox',d,1),'edgecolor','interp')
    colorbar
    caxis([floor(min(min(F_SDEtoolbox(501:end,:)))) floor(max(max(F_SDEtoolbox(501:end,:))))])
    caxis([0 2])
    colormap(gray)
    xlim([0 tspan(end)-T_initial])
    xlabel('time','fontsize',20)
    ylabel('cell perimeter','fontsize',20)
    title(['Figure 7F: F Kymograph'],'fontsize',20,'fontweight','n')
%%
% R kymograph 
    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    surf(T-T_initial,x,circshift(R_SDEtoolbox',d,1),'edgecolor','interp')
    colorbar

caxis([3 8])
    colormap(gray)
    xlim([0 tspan(end)-T_initial])
    xlabel('time','fontsize',20)
    ylabel('cell perimeter','fontsize',20)
    title('Figure 7F: R Kymograph','fontsize',20,'fontweight','n')
%%
% B kymograph 

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
    title('Figure 7F: B Kymograph','fontsize',20,'fontweight','n')


%%
%   Cpm kymograph
    Cpm_SDEtoolbox = Cpmf_SDEtoolbox + Cpmb_SDEtoolbox;
    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    surf(T-T_initial,x,circshift(Cpm_SDEtoolbox',d,1),'edgecolor','interp')
    colorbar
    colormap(gray)
    xlim([0 tspan(end)-T_initial])
    xlabel('time','fontsize',20)
    ylabel('cell perimeter','fontsize',20)
    title('Figure 7F: $$C^+_m =  C^+_{mf}+ C^+_{mb} \quad \mathrm{Kymograph}$$ ','Interpreter','Latex','fontweight','n','fontsize',20)
    caxis([20 50])
    
%%

    sp_index = 200;
    v1 = F_SDEtoolbox(:,sp_index);
    v2 = B_SDEtoolbox(:,sp_index);
    v3 = Cpmb_SDEtoolbox(:,sp_index);
    v4 = Cpm_SDEtoolbox(:,sp_index);

    v1 = v1(1:10:end);
    v2 = v2(1:10:end);
    v3 = v3(1:10:end);
    v4 = v4(1:10:end);

    t1 = T-T_initial;
    t1 = t1(1:10:end);


    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    plot(t1,smooth(v1/max(v1)),'g','linewidth',2)
    plot(t1,smooth(v2/max(v2)),'r','linewidth',2)
    plot(t1,smooth(v4/max(v4)),'-','color',[255,140,0]/255,'linewidth',2)
    xlim([0 600])
    xlabel('Time (A.U.)')
    ylabel('Normalized Intensity')
    title(['Figure 7G Top'],'fontweight','n','fontsize',20)
    
%%

    sp_index = 240;
    v1 = F_SDEtoolbox(:,sp_index);
    v2 = B_SDEtoolbox(:,sp_index);
    v3 = Cpmb_SDEtoolbox(:,sp_index);
    v4 = Cpm_SDEtoolbox(:,sp_index);

    v1 = v1(1:10:end);
    v2 = v2(1:10:end);
    v3 = v3(1:10:end);
    v4 = v4(1:10:end);

    t1 = T-T_initial;
    t1 = t1(1:10:end);


    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    plot(t1,smooth(v1/max(v1)),'g','linewidth',2)
    plot(t1,smooth(v2/max(v2)),'r','linewidth',2)
    plot(t1,smooth(v4/max(v4)),'-','color',[255,140,0]/255,'linewidth',2)
    xlim([0 600])
    xlabel('Time (A.U.)')
    ylabel('Normalized Intensity')
    title(['Figure 7G Bottom'],'fontweight','n','fontsize',20)
    

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
    


    F_val = linspace(0.01,6,10001);
    B_val_0 = param.b1./(param.b3.*F_val + 0*param.Inh*10);
    B_val_1 = param.b1./(param.b3.*F_val + param.Inh*10);
    
    yR = (param.c2.*F_val.*param.C2)./(param.C1.*param.c1+param.c2.*F_val);
    yF_p_0 = (-param.a1.*F_val + ((param.A1.*param.a3+ param.phi*3.86)./(param.A2.*param.a4^2.*B_val_0.^2+1)+param.ub).*(param.a5-F_val))./(param.a2.*F_val);
    yF_np_0 = (-param.a1.*F_val + ((param.A1.*param.a3+ param.phi*-1.2)./(param.A2.*param.a4^2.*B_val_0.^2+1)+param.ub).*(param.a5-F_val))./(param.a2.*F_val);
    
    yF_np_1 = (-param.a1.*F_val + ((param.A1.*param.a3+ param.phi*-1.2)./(param.A2.*param.a4^2.*B_val_1.^2+1)+param.ub).*(param.a5-F_val))./(param.a2.*F_val);
    

    
    
    figure;
    hold on
    plot(F_val,B_val_0)
    plot(F_val,B_val_1)
    plot(F_val,yR)

    plot(F_val,yF_np_0)
    plot(F_val,yF_np_1)
    xlim([0 3])
    ylim([0 10])
end

