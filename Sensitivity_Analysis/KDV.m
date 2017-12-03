% @author: Maziar Raissi

function params_list = KDV()
% quantile(params_list,[0.025 0.25 0.50 0.75 0.975])
clc; close all;

plt = 1;
save_plt = 0;

addpath ..
addpath ../Utilities
addpath ../Kernels/KDV
addpath ../Utilities/export_fig

function CleanupFun()
    rmpath ..
    rmpath ../Utilities
    rmpath ../Kernels/KDV
    rmpath ../Utilities/export_fig
end

finishup = onCleanup(@() CleanupFun());

rng('default')

set(0,'defaulttextinterpreter','latex')

%% Load Data
load('../Data/kdv.mat', 'usol', 't', 'x')
u_star = real(usol); % 512x201
t_star = t; % 201x1
x_star = x';   % 512x1
N_star = size(x_star,1);
nsteps = size(t_star,1)-1;

if plt ==1
    figure(1);
    plot_surface(t_star, x_star, u_star, '$t$', '$x$', '$u(t,x)$');
    view(3)

    drawnow()
end
    
%% Setup
noise = 0.00;
u_data = u_star + noise*std(u_star(:))*randn(size(u_star));

N0 = 111;
N1 = 109;
%% Optimize model
params_list = zeros(nsteps,2);
hyp = [log([1.0 1.0]) 0.0 0.0 -4.0];
idx1 = randsample(N_star, N0);
step = 1;
for i = 1:step:nsteps
    dt = t_star(i+step) - t_star(i);
    
    idx0 = idx1;
    x0 = x_star(idx0,:);
    u0 = u_data(idx0,i);
    
    idx1 = randsample(N_star,N1);
    x1 = x_star(idx1,:);
    u1 = u_data(idx1,i+step);
    
    model = HPM(x1, u1, x0, u0, dt, hyp);
    model = model.train(50);
    
    hyp = model.hyp;
    params_list(i,:) = hyp(3:4);
    
    [pred_n_star, var_n_star] = model.predict(x_star);
    var_n_star = abs(diag(var_n_star));
    
    error = norm(pred_n_star - u_star(:,i+step))/norm(u_star(:,i+step));
    
    fprintf(1,'=========================\n');
    fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error = %.2e\n\n', i, ...
        t_star(i+step), model.NLML, error);
       
    str = sprintf('%.2f  ', params_list(i,:));
    fprintf('Parameters: %s\n\n', str)
    
    str = sprintf('%.2f  ', median(params_list(1:step:i,:),1));
    fprintf('Median: %s\n', str)
    fprintf(1,'=========================\n\n');
    
    if plt == 1
        if ~exist('fig','var')
            fig = figure(2);
        end
        set(fig,'units','normalized','outerposition',[0 0 1 1])
        clf
        
        subplot(3,1,1);
        tit = sprintf('Time: %.2f\n%d training points', t_star(i), N0);
        plot_data_1D(x_star, u_star(:,i), x0, u0, '$x$', '$u(t,x)$', tit);
        
        
        subplot(3,1,2);
        tit = sprintf('Time: %.2f\n%d training points', t_star(i+step), N1);
        plot_data_1D(x_star, u_star(:,i+step), x1, u1, '$x$', '$u(t,x)$', tit);
        
        
        subplot(3,1,3);
        plot_prediction_1D(x_star, u_star(:,i+step), pred_n_star, var_n_star, ...
            '$x$', '$u(t,x)$', tit);
        
        drawnow;
    end    
    
end

if save_plt == 1
    export_fig ./Figures/KDV.png -r300
end

end