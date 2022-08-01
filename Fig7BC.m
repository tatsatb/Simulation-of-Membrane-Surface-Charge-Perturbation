function Fig7BC

    clear 
    close all
    clc
    warning('off','all');  

    global param

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

findeqm

end

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
    

    
    
    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    plot(F_val,B_val_0,'r-','LineWidth',2)
    plot(F_val,B_val_1,'r--','LineWidth',2)
    legend('without actuation','with actuation','fontsize',20)
    xlim([0 3])
    ylim([0 10])
    xlabel('F')
    ylabel('B')
    title('Figure 7B','FontWeight','normal')


    figure('color','white')
    set(gca, 'fontweight','n','linewidth',1,'fontsize',20)
    hold on
    plot(F_val,yR,'b-','LineWidth',2)
    h(1) = plot(F_val,yF_np_0,'g-','LineWidth',2);
    h(2) = plot(F_val,yF_np_1,'g--','LineWidth',2);
    legend(h,{'without actuation','with actuation'},'fontsize',20)
    xlim([0 3])
    ylim([0 10])
     xlabel('F')
    ylabel('R')
    title('Figure 7C','FontWeight','normal')
end