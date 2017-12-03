% @author: Maziar Raissi

function params_list = Navier_Stokes()
% quantile(params_list,[0.025 0.25 0.50 0.75 0.975])
clc; close all;

plt = 1;
save_plt = 0;

addpath ..
addpath ../Utilities
addpath ../Kernels/Navier_Stokes
addpath ../Utilities/export_fig

function CleanupFun()
    rmpath ../Utilities
    rmpath ../Kernels/Navier_Stokes
    rmpath ../Utilities/export_fig
end

finishup = onCleanup(@() CleanupFun());

rng('default')

set(0,'defaulttextinterpreter','latex')


%% Load Data
load('../Data/cylinder_fine.mat', 't_star', 'X_star', 'U_star', 'w_star')
N_star = size(X_star,1);
nsteps = size(t_star,1) - 1;

    
%% Setup
noise = 0.00;
U_data = U_star;
u_star = U_star(:,1,:);
U_data(:,1,:) = U_star(:,1,:) + noise*std(u_star(:))*randn(size(U_star(:,1,:)));
v_star = U_star(:,2,:);
U_data(:,2,:) = U_star(:,2,:) + noise*std(v_star(:))*randn(size(U_star(:,2,:)));

N0 = 251;
N1 = 249;
%% Optimize model
params_list = zeros(nsteps,2);
hyp = [log([1.0 1.0 1.0]) 0.0 -1.0 log([1.0 1.0 1.0]) -4.0];
idx1 = randsample(N_star, N0);
step = 1;
for i = 1:step:nsteps
    dt = t_star(i+step) - t_star(i);
    
    idx0 = idx1;
    X0 = X_star(idx0,:);
    U0 = U_data(idx0,:,i);
    
    idx1 = randsample(N_star,N1);
    X1 = X_star(idx1,:);
    U1 = U_data(idx1,:,i+step);
    
    model = HPM(X1, U1, X0, U0, dt, hyp);
    model = model.train(20);
    
    hyp = model.hyp;
    params_list(i,:) = [hyp(4) exp(hyp(5))];
        
    [pred_n_star, var_n_star] = model.predict(X_star);
    pred_U_star = reshape(pred_n_star,N_star,2);
    var_n_star = abs(diag(var_n_star));
    
    error = norm(pred_U_star(:,1) - U_star(:,1,i+step))/norm(U_star(:,1,i+step));
    
    fprintf(1,'=========================\n');
    fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error = %.2e\n\n', i, ...
        t_star(i+step), model.NLML, error);
       
    str = sprintf('%.4f  ', params_list(i,:));
    fprintf('Parameters: %s\n\n', str)
    
    str = sprintf('%.4f  ', median(params_list(1:step:i,:),1));
    fprintf('Median: %s\n', str)
    fprintf(1,'=========================\n\n');
    
    if plt == 1
        if ~exist('fig','var')
            fig = figure(2);
        end
        set(fig,'units','normalized','outerposition',[0 0 1 1])
        clf
        
        subplot(2,2,1)
        plot_surface_griddata(X_star, U_star(:,1,i+step),'','','')
        
        subplot(2,2,2);
        plot_surface_griddata(X_star, U_star(:,2,i+step),'','','')
        
        subplot(2,2,3)
        plot_surface_griddata(X_star, pred_U_star(:,1),'','','')
        
        subplot(2,2,4);
        plot_surface_griddata(X_star, pred_U_star(:,2),'','','')       

        drawnow()
    end    
    
end

if save_plt == 1
    export_fig ./Figures/Navier_Stokes.png -r300
end


end