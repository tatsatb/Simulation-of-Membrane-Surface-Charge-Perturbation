

function [out1,out2,out3] = Fig7FG_sdefile(t,x,flag,initvalues,SDETYPE,NUMDEPVARS,NUMSIM,param,input)

    global param

% Initial condition
    Fzero       = initvalues(1);     % F
    Rzero       = initvalues(2);     % R
    Bzero       = initvalues(3);     % B

    Wzero       = initvalues(4);     
    Zzero       = initvalues(5);     

    
    Cpmfzero    = initvalues(6);     
    Cpmbzero    = initvalues(7);   

% IN = input.baseMean*ones(1,param.Np);

if nargin < 3 || isempty(flag)
    
    xsplitted  =  SDE_split_sdeinput(x,NUMDEPVARS);
   
  %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
  %:::::::::::::::::::::::::::::::  DEFINE HERE THE SDE  :::::::::::::::::::::::::::::
  %::::::::::::: (define the initial conditions at the bottom of the page) :::::::::::
  %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    F         = xsplitted{1}; % F
    R         = xsplitted{2}; % R
    B         = xsplitted{3}; % R
  
    W         = xsplitted{4}; 
    Z         = xsplitted{5}; 
  
 
    Cpmf      = xsplitted{6}; 
    Cpmb      = xsplitted{7}; 
  
  

           
   switch upper(SDETYPE)
   case 'ITO'
      

   
   if t>30
       P_switch = 1;
   else
       P_switch = 0;
   end
   
   
   if t > 180 && t<= 360
            
           Cp_switch1 = 0*F;
           Cp_switch1(1:end) = 1; % global
%            Cp_switch1(40:90) = 1; % local
                
   else
           Cp_switch1 = 0;
   end  


  

%=======================================================       

       
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        Fshiftleft  =  [F(:,end) F(:,1:end-1)];
        Fshiftright =  [F(:,2:end) F(:,1)];
        

         driftF      =    param.T_scale*( -(param.a1+param.a2*R).*F + ((param.A1.*param.a3+ param.phi*(W-Z).*P_switch)./(param.A2.*param.a4^2.*B.^2+1)+param.ub).*(param.a5-F)) ...
                       +  param.D_F.*(Fshiftleft + Fshiftright - 2*F);
         diffusionF  =    1*sqrt(abs(param.T_scale*( -(param.a1+param.a2*R).*F + ((param.A1.*param.a3)./(param.A2.*param.a4^2.*B.^2+1)+param.ub).*(param.a5-F)))) ...
                       +  0*input.baseSigma;           
        derivativeF = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Rshiftleft  =  [R(:,end) R(:,1:end-1)];
        Rshiftright =  [R(:,2:end) R(:,1)];
             
        driftR      =   param.T_scale*(-param.C1.*param.c1.*R + param.c2.*F.*(param.C2-R))...
                      + param.D_R.*(Rshiftleft + Rshiftright - 2*R);  
  
        diffusionR =  1 *sqrt(abs(param.T_scale*(-param.C1*param.c1*R + param.c2*F.*4)));
        derivativeR = 0; 
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
        Bshiftleft  =  [B(:,end) B(:,1:end-1)];
        Bshiftright =  [B(:,2:end) B(:,1)];

        driftB      =   param.T_scale*100*(param.b1 - param.b3*F.*B - param.Inh.*Cpmb.*B)...
                      + param.D_B.*(Bshiftleft + Bshiftright - 2*B);
        diffusionB =  1*sqrt(abs(100*(param.T_scale*(param.b1 - param.b3*F.*B - param.Inh.*Cpmb.*B))));
        derivativeB = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
        Wshiftleft  =  [W(:,end)   W(:,1:end-1)];
        Wshiftright =  [W(:,2:end) W(:,1)];

        driftW      =   param.T_scale*(5*(R-W)) ...
                         + 1*param.D_W.*(Wshiftleft + Wshiftright - 2*W);
        diffusionW =  sqrt(abs(param.T_scale*(5*(R-W))));
        derivativeW = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
        Zshiftleft  =  [Z(:,end)   Z(:,1:end-1)];
        Zshiftright =  [Z(:,2:end) Z(:,1)];

%     
        driftZ      =   param.T_scale*(1*(R-Z)) ...
                         + 1*param.D_Z.*(Zshiftleft + Zshiftright - 2*Z);             
        diffusionZ =  sqrt(abs(param.T_scale*(1*(R-Z))));
        derivativeZ = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
      
        Cpmfshiftleft  =  [Cpmf(:,end)   Cpmf(:,1:end-1)];
        Cpmfshiftright =  [Cpmf(:,2:end) Cpmf(:,1)];
                     
        driftCpmf      =   param.T1*param.T_scale*(+Cp_switch1*100 - 0.1.*Cpmf.*B + 8.*Cpmb.*F-3*Cpmf)...
                         + param.D_Cpmf.*(Cpmfshiftleft + Cpmfshiftright - 2*Cpmf);
        diffusionCpmf =  1*sqrt(abs(param.T1*param.T_scale*(+Cp_switch1*100 - 0.1.*Cpmf.*B + 8.*Cpmb.*F-3*Cpmf)));
        derivativeCpmf = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
        Cpmbshiftleft  =  [Cpmb(:,end)   Cpmb(:,1:end-1)];
        Cpmbshiftright =  [Cpmf(:,2:end) Cpmb(:,1)];

                     
        driftCpmb      =   param.T2*param.T_scale*(0.1.*Cpmf.*B - 8.*Cpmb.*F) ...
                         + 0*param.D_Cpmb.*(Cpmbshiftleft + Cpmbshiftright - 2*Cpmb);
        diffusionCpmb =  1*sqrt(abs(param.T2*param.T_scale*(0.1.*Cpmf.*B - 8.*Cpmb.*F)));
        derivativeCpmb = 0; 
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
     
        
     
        
   end
   
    out1 = zeros(1,NUMDEPVARS*NUMSIM);
    out1(1:NUMDEPVARS:end) = driftF;
    out1(2:NUMDEPVARS:end) = driftR;
    out1(3:NUMDEPVARS:end) = driftB;
  
    out1(4:NUMDEPVARS:end) = driftW;
    out1(5:NUMDEPVARS:end) = driftZ;
  

    out1(6:NUMDEPVARS:end) = driftCpmf;
    out1(7:NUMDEPVARS:end) = driftCpmb;
%--------------------------------------------------------------------------    
    out2 = zeros(1,NUMDEPVARS*NUMSIM);
    out2(1:NUMDEPVARS:end) = diffusionF;
    out2(2:NUMDEPVARS:end) = diffusionR;
    out2(3:NUMDEPVARS:end) = diffusionB;
   
    out2(4:NUMDEPVARS:end) = diffusionW;
    out2(5:NUMDEPVARS:end) = diffusionZ;

    out2(6:NUMDEPVARS:end) = diffusionCpmf;
    out2(7:NUMDEPVARS:end) = diffusionCpmb;
%--------------------------------------------------------------------------        
    out3 = zeros(1,NUMDEPVARS*NUMSIM);
    out3(1:NUMDEPVARS:end) = derivativeF;
    out3(2:NUMDEPVARS:end) = derivativeR;
    out3(3:NUMDEPVARS:end) = derivativeB;
   
    out3(4:NUMDEPVARS:end) = derivativeW;
    out3(5:NUMDEPVARS:end) = derivativeZ;

    out3(6:NUMDEPVARS:end) = derivativeCpmf;
    out3(7:NUMDEPVARS:end) = derivativeCpmb;
    
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   %:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
else
    
    switch(flag)
    case 'init'  
        out1 = t;
        
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%::::::::::::::::::::::  DEFINE HERE THE SDE INITAL CONDITIONS  :::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

out2 = [Fzero Rzero Bzero Wzero Zzero Cpmfzero Cpmbzero];   % write here the SDE initial condition(s)
        
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: :::::::::::::::::::::
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        out3 = [];
        
        
    otherwise
        error(['Unknown flag ''' flag '''.']);
    end
end
