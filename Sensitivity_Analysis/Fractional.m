% @author: Maziar Raissi

function params_list = Fractional()
% quantile(params_list,[0.025 0.25 0.50 0.75 0.975])

clc; close all;

plt = 1;
save_plt = 0;

addpath ..
addpath ../Utilities
addpath ../Kernels/Fractional
addpath ../Utilities/export_fig

function CleanupFun()
    rmpath ..
    rmpath ../Utilities
    rmpath ../Kernels/Fractional
    rmpath ../Utilities/export_fig
end

finishup = onCleanup(@() CleanupFun());

rng('default')

set(0,'defaulttextinterpreter','latex')

%% Load Data
[t_star, x_star, u_star] = gen_data(10^6, 0.01, 5, 100);
% u_star --> 300x5
% t_star --> 5x1
% x_star --> 300x1
N_star = size(x_star,1);
nsteps = size(t_star,1)-1;

if plt ==1
    figure(2);
    plot_surface(t_star, x_star, u_star, '$t$', '$x$', '$u(t,x)$');
    view(3);

    drawnow()
end
    
%% Setup
noise = 0.0;
u_data = (1 + noise*randn(size(u_star))).*u_star;

N0 = 100;
N1 = 100;
%% Optimize model
params_list = [];
hyp = [log([1.0 0.1]) 0.1 1.2 -4.0];
idx1 = randsample(N_star, N0);
step = 1;
for i = nsteps:step:nsteps
    dt = t_star(i+step) - t_star(i);
    
    idx0 = idx1;
    x0 = x_star(idx0,:);
    u0 = u_data(idx0,i);
    
    idx1 = randsample(N_star,N1);
    x1 = x_star(idx1,:);
    u1 = u_data(idx1,i+step);
    
    model = HPM(x1, u1, x0, u0, dt, hyp);
    model = model.train(500);
    
    hyp = model.hyp;
    params_list = [params_list; hyp(3:4)];
    
    [pred_n_star, var_n_star] = model.predict(x_star);
    var_n_star = abs(diag(var_n_star));
    
    error = norm(pred_n_star - u_star(:,i+step))/norm(u_star(:,i+step));
    
    fprintf(1,'=========================\n');
    fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error = %.2e\n\n', i, ...
        t_star(i+step), model.NLML, error);
       
    str = sprintf('%.2f  ', params_list(end,:));
    fprintf('Parameters: %s\n\n', str)
    
    str = sprintf('%.2f  ', median(params_list,1));
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
    export_fig ./Figures/Burgers.png -r300
end

end

function [t_star, x_star, u_star] = gen_data(N,dt,m,n)

    pos = cumsum(sqrt(dt)*randn(N,1));

    M = 0;

    P = zeros(length(pos)-m,m);
    for i = 1:length(pos)-m
        y = pos(i+1:i+m) - pos(i);
        M = max([M, max(abs(y))]);
        P(i,:) = y;
    end
    
    M = M/2;
    
    bins = linspace(-M,M,n+1);
    x_star = linspace(M*(1/n-1), M*(1-1/n), n)';
    t_star = dt*(1:1:m)';
    u_star = zeros(n,m);
    
    figure(1)
    clf
    hold
    for i = 1:m
        h = histogram(P(:,i), bins, 'Normalization', 'pdf');
        u_star(:,i) = h.Values;
    end

end
