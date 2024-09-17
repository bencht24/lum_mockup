
% NOTE: Script takes significant influence form Lehmann et al.'s OSL MLE script,
% accessible on github: https://github.com/BenjaminLehmann/Esurf2019_new_codes
% differences include using LAR protocol, the removal of calculating confidence intervals, removing log
% likelihood in plot (instead plotting inverse error for relative
% view), and plotting the PDF for reference of shape only
clear all; 
clc;
close all;

TT  = 10000;

% Loading data file excel
SampleName     = 'OSL_Depth_Data_agecalc.xlsx';
[num]          = xlsread(SampleName);
n              = length(num(:,1));

% Input parameters from prior function
SP0 = 3.20E-07;                      % [s-1]
mu  = .285;                        % [mm-1]
SP  = SP0*365.25*24.*3600;           % [a-1]

ind_plt      =24;                   % row # where depth profile reaches saturation (curve plateau)

%% Parameterization of time
% Enter time range to be sampled. 

tmin = 0;
tmax = 1000;

% Definitions de donnees experimentales

x(1:n)          = (num(:,1)); % first column from excel sheet, depth [mm]
L(1:n)          = (num(:,2));  % second column from excel sheet, Lx/Tx
e(1:n)          = (num(:,3)); % third column from excel sheet, error of Lx/Tx
[x_s,ind]       = sort(x(:)); % sorting increasingly the depth, if not done already
Ls_M              = L(ind); % sorting Lx/Tx with increasing depth
a               = std(Ls_M(ind_plt:n)); % calculating the standard deviation on the plateau

%% Compute residuals (i.e. fit to data)

M         = nan(TT,1);
t_vec     = nan(TT,1);

h         = waitbar(0,'Thanks Benjamin Lehmann for publishing the script I based this protocol on!');  % creation of the loading bar (just for fun)
rand_vec  = rand(TT,1);

r_t1      = sort(tmin+(tmax-tmin)*rand_vec);

for i = 1:TT
            
            M(i)   = 0;           
            L_th   = exp(-SP*r_t1(i)*exp(-mu*x));        % Equation for misfit calc generated
            M(i)   = sum(((L-L_th).^2)./a);             % Calculate the misfit M          
            t_vec(i)  = r_t1(i);

        waitbar(i/TT,h)
end
    
close(h)

chi      = 1./exp(0.5*M);    % Likelihood non normalized
max_chi  = max(chi(:));      % Max value of the Likelihood
norm_chi = chi/max_chi;      % Likelihood normalized

%% (Archaic) was used originally to resampling the likelihood by comparing to a random value between 0 and 1

jt=0;
for it=1:TT
    R=rand;
    if (M(it)>-1)
        jt        = jt+1;
        s_chi(jt) = M(it);
        s_t(jt)   = t_vec(it);
    end
end
% old sequence, kept in for consistency

%locates best fit from matrix
 [minValue, linearIndexesOfMaxes] = min(M(:));
[rowsOfMaxes colsOfMaxes] = find(M == minValue);   
   mode_bestfit = M(rowsOfMaxes,colsOfMaxes);
    r_t_bestfit = t_vec(rowsOfMaxes);
    disp(['best fit - age = ' num2str(r_t_bestfit,10)]);



%% extract misfit
LxTx_BF  = exp(-SP0*r_t_bestfit*exp(-mu*x));    % creation of synthetic luminescence signal
             %ERR   = (((L-L_th).^2)./a); 
            ERR   = (abs(L-L_th))./a;   % calculation of the misfit between experimental data and synthetic signal
           % ((L-L_th).^2)./a
                  %writematrix(ERR,'misfitage.txt'); %best fit misfit error export of each depth profile datum


%% Plotting the distribution to showcase viability of best fit (reference)

  
plot(t_vec,norm_chi,'b','LineWidth',1)
hold on

