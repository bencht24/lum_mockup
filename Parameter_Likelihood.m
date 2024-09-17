
% NOTE: Script takes significant influence form Lehmann et al.'s OSL MLE script,
% accessible on github: https://github.com/BenjaminLehmann/Esurf2019_new_codes
% differences include the removal of calculating confidence intervals, removing log
% likelihood in surface plot (instead plotting inverse error for relative
% view), and plotting best fit point to surface plot data. 
 
clear all;
clc;
close all;


TT      = 1000;       % Size of inversion matrix used for tested parameter

%% INput data from file excel 

[num,txt,tab]    =  xlsread('OSL_Depth_Data.xlsx'); %be sure to include data from all cores!
n                =  length(num(:,1));

%% Input of the known age in years.

t_known          = 11;    %0.031709792  % Input the known age of the sample in years
t                = t_known*365.25*24.0*3600.0;  % Conversion in seconds

%% Ordering of the imported data

Coredat          = 32;                     % Number of data --> last row # of data on excel sheet
sat_pt      =24;                     % row # where depth profile reaches saturation (curve plateau)

x(1:n)       = (num(1:n,1));           % first column from excel sheet, depth [mm]
LxTx(1:n)    = (num(1:n,2));           % second column from excel sheet, Lx/Tx
LxTx_e(1:n)  = (num(1:n,3));           % third column from excel sheet, error of Lx/Tx

[x_s,ind]    = sort(x(:));             % sorting increasingly the depth, if not done already
LxTx_s       = LxTx(ind);              % sorting Lx/Tx with increasing depth
LxTx_a       = std(LxTx_s(sat_pt:n)); % calculating the standard deviation on the plateau

%% Define tested parameter range 


log_SP_t_max = 10; 
log_SP_t_min = 0;     
SP_t_max     = 10^(log_SP_t_max);   % Dimensionless
SP_t_min     = 10^(log_SP_t_min);   % Dimensionless
         
mu_max       = 1;                 % mm-1
mu_min       = 0;                 % mm-1

%% Initialization of the variable and matrix for the inversion

M            = nan(TT,TT);    % creation of an empty misfit matrix

mu_matrix    = nan(TT,TT);    % creation of an empty mu matrix
SP_t_matrix  = nan(TT,TT);    % creation of an empty SP_t matrix

rand_vec1    = rand(TT,1);    % random vector of size TT to sample SP_t
rand_vec2    = rand(TT,1);    % random vector of size TT to sample mu

r_SP_t       = sort(10.^(log_SP_t_min+(log_SP_t_max-log_SP_t_min)*rand_vec1));    % r_SP_t is randomly sampled in log space
r_mu         = sort(mu_min+(mu_max-mu_min)*rand_vec2);                            % mu is randomly sampled in normal space

%% Inversion

h               = waitbar(0,'Thanks to Benjamin Lehmann for publishing the script I based this protocol on!');  % creation of the loading bar (just for fun)

    for i = 1:TT
        for j = 1:TT                                       % "parfor" for paralellelization, for normal claculation put "for"
            
            M(i,j)   = 0;                                  % initialization of the misfit at i and j index
                                
            LxTx_th  = exp(-r_SP_t(j)*exp(-r_mu(i)*x));    % creation of synthetic luminescence signal
            M(i,j)   = sum((abs(LxTx-LxTx_th))./LxTx_a);   % calculation of the misfit between experimental data and synthetic signal                             
                                    
            mu_matrix(i,j)  = r_mu(i);                     % compile mu tested into matrix
            SP_t_matrix(i,j)  = r_SP_t(j);                 % compile SP_t tested into matrix
            
        end
        waitbar(i/TT,h)                                    % incrementation of the loading bar (just for fun)
    end

close(h)                                                   % close the loading bar

%% Transformation of the misfit M into likelihood chi

chi_m        = 1./exp(0.5*M);

%% (archaic) Originally used for resampling the likelihood by comparing to a random value between 0 and 1
%made to incorporate all fits into plot, --> M(it,j) >-1. This is because
%PDF's are not used for certainty in this script. Rather an inverse jacobian function is used, see A3 of appendix 
jt=0; 
for it=1:TT
    for j= 1:TT
    R=rand;
    if (M(it,j)>-1)
        jt=jt+1;
        s_chi(jt)= M(it,j); 
        s_SP_t(jt) = SP_t_matrix(it,j);
        s_mu(jt) = mu_matrix(it,j);
    end
    end 
end


chi        = 1./exp(0.5*s_chi); %used to identify likelihood best fit parameters from matrix, which is why this is above resampling code is kept in the script




%% Print the results

fprintf('\nResult for the calibration\n')


%% Creating output models

%LxTx_M   = exp(-SP_t_M_ok*exp(-mu_M*x_model));

%% Convert matrix of SP_t [dimensioneless] into sigmaphio parameter SP [s-1]

SP_max     = SP_t_max/t;   % s-1
SP_min     = SP_t_min/t;   % s-1

r_SP       = r_SP_t./t;    % s-1 
%% Identify best fitting mu and sigmaphio parameters from matrix
 [minValue, linearIndexesOfMaxes] = min(M(:));
[rowsOfMaxes colsOfMaxes] = find(M == minValue);   
   mode_bestfit = M(rowsOfMaxes,colsOfMaxes);
    r_mu_bestfit = r_mu(rowsOfMaxes);
    r_sig_bestfit = r_SP(colsOfMaxes);
    disp(['best fit - attenuation = ' num2str(r_mu_bestfit,5)]);
disp(['best fit - bleaching rate = ' num2str(r_sig_bestfit,5)]);
r= max(chi);

%%produce fit error quantity for best fit for each datum
LxTx_BF  = exp(-r_sig_bestfit*t*exp(-r_mu_bestfit*x));    % creation of synthetic luminescence signal
            ERR   = (abs(LxTx-LxTx_th))./LxTx_a;   % calculation of the misfit between experimental data and synthetic signal                             
                  writematrix(ERR,'misfit.txt')   %error provided for each datum using best fit parameters
%% Generate surface plot

surface(r_SP,r_mu,chi_m); axis square; colorbar;shading interp; %surface plots sigphio, mu, and resultant inverse misfit
set(gca,'XScale','log');
title(colorbar,'Likelihood (Inverse Misfit)')
hold on
p = scatter3(r_sig_bestfit,r_mu_bestfit,r, 10, 'filled', 'r'); %plots best fit point

hold off

xlabel('Bleaching rate [s-1]')
ylabel('Attenuation coeff. [mm-1]')
axis([SP_min SP_max mu_min mu_max])
legend([p],{'Best Fit'});
title('Probability distribution')