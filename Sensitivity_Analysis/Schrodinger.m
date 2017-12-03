% @author: Maziar Raissi

function params_list = Schrodinger()
% quantile(params_list,[0.025 0.25 0.50 0.75 0.975])
clc; close all;

plt = 1;
save_plt = 0;

addpath ..
addpath ../Utilities
addpath ../Kernels/Schrodinger
addpath ../Utilities/export_fig

function CleanupFun()
    rmpath ..
    rmpath ../Utilities
    rmpath ../Kernels/Schrodinger
    rmpath ../Utilities/export_fig
end

finishup = onCleanup(@() CleanupFun());

rng('default')

set(0,'defaulttextinterpreter','latex')

%% Load Data
load('../Data/nls.mat', 'usol', 't', 'x')
u_star = real(usol)'; % 501x512
v_star = imag(usol)'; % 501x512
t_star = t; % 501x1
x_star = x';   % 512x1
N_star = size(x_star,1);
nsteps = size(t_star,1)-1;

if plt ==1
    fig = figure(1);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    subplot(3,1,1);
    plot_surface(t_star, x_star, u_star, '$t$', '$x$', '$u(t,x)$');
    view(3)
    
    subplot(3,1,2);
    plot_surface(t_star, x_star, v_star, '$t$', '$x$', '$v(t,x)$');
    view(3)
    
    subplot(3,1,3);
    plot_surface(t_star, x_star, (u_star.^2 + v_star.^2).^(1/2), '$t$', '$x$', '$(u(t,x)^2 + v(t,x)^2)^{1/2}$');
    view(3)
    
    drawnow()
end
    
%% Setup
noise = 0.00;
u_data = u_star + noise*std(u_star(:))*randn(size(u_star));
v_data = v_star + noise*std(v_star(:))*randn(size(v_star));

N0 = 49;
N1 = 51;
%% Optimize model
params_list = zeros(nsteps,2);
hyp = [log([1.0 1.0]) log([1.0 1.0]) 0.0 0.0 -4.0];
idx1 = randsample(N_star, N0);
step = 1;
for i = 1:step:nsteps
    dt = t_star(i+step) - t_star(i);
    
    idx0 = idx1;
    x0 = x_star(idx0,:);
    U0 = [u_data(idx0,i) v_data(idx0,i)];
    
    idx1 = randsample(N_star,N1);
    x1 = x_star(idx1,:);
    U1 = [u_data(idx1,i+step) v_data(idx1,i+step)];
    
    model = HPM(x1, U1, x0, U0, dt, hyp);
    
    
    model = model.train(50);
    
    
    hyp = model.hyp;
    params_list(i,:) = hyp(5:6);
    
    [pred_n_star, var_n_star] = model.predict(x_star);
    pred_U_star = reshape(pred_n_star,N_star,2);
    var_n_star = abs(diag(var_n_star));
    var_U_star = reshape(var_n_star,N_star,2);
    
    error_u = norm(pred_U_star(:,1) - u_star(:,i+step))/norm(u_star(:,i+step));
    error_v = norm(pred_U_star(:,2) - v_star(:,i+step))/norm(v_star(:,i+step));
    
    fprintf(1,'=========================\n');
    fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error_u = %.2e, Error_v = %.2e\n\n', i, ...
        t_star(i+step), model.NLML, error_u, error_v);
       
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
        plot_data_1D(x_star, u_star(:,i), x0, U0(:,1), '$x$', '$u(t,x)$', tit);
        
        
        subplot(3,1,2);
        tit = sprintf('Time: %.2f\n%d training points', t_star(i+step), N1);
        plot_data_1D(x_star, u_star(:,i+step), x1, U1(:,1), '$x$', '$u(t,x)$', tit);
        
        
        subplot(3,1,3);
        plot_prediction_1D(x_star, u_star(:,i+step), pred_U_star(:,1), var_U_star(:,1), ...
            '$x$', '$u(t,x)$', tit);
        
        drawnow;
    end    
    
end

if save_plt == 1
    export_fig ./Figures/Burgers.png -r300
end

end