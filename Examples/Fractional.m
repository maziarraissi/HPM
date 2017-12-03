% @author: Maziar Raissi

function params = Fractional()

clc; close all;

plt = 1;
plt_pred = 0;
save_plt = 1;

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
[t_star, x_star, u_star, pos] = gen_data(10^6, 0.01, 5, 100);
% u_star --> 100x5
% t_star --> 5x1
% x_star --> 100x1
N_star = size(x_star,1);
nsteps = size(t_star,1)-1;
    
%% Setup
N0 = 100;
N1 = 100;
%% Clean Data
i = 4;
dt = t_star(i+1) - t_star(i);

idx0 = randsample(N_star, N0);
x0 = x_star(idx0,:);
u0 = u_star(idx0,i);

idx1 = randsample(N_star,N1);
x1 = x_star(idx1,:);
u1 = u_star(idx1,i+1);
    
hyp = [log([1.0 0.1]) 0.1 1.2 -4.0];
model = HPM(x1, u1, x0, u0, dt, hyp);
model = model.train(500);
    
hyp = model.hyp;
params = hyp(3:4);
    
[pred_n_star, var_n_star] = model.predict(x_star);
var_n_star = abs(diag(var_n_star));
    
error = norm(pred_n_star - u_star(:,i+1))/norm(u_star(:,i+1));
    
fprintf(1,'=========================\n');
fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error = %.2e\n\n', i, ...
    t_star(i+1), model.NLML, error);

str = sprintf('%.2f  ', params);
fprintf('Parameters: %s\n\n', str)
fprintf(1,'=========================\n\n');
    
if plt_pred == 1
    figure();
    plot_prediction_1D(x_star, u_star(:,i+1), pred_n_star, var_n_star, ...
        '$x$', '$u(t,x)$', 'Prediction');
    
    drawnow;
end

%% Plot Results

if plt == 1
    fig = figure();
    set(fig,'units','normalized','outerposition',[0 0 1 .5])
    subplot(3,2,1:2)
    plot(pos, 'b', 'LineWidth', 2);
    xlabel('Time')
    ylabel('Position')
    axis tight
    set(gca,'FontSize',14);
    set(gcf, 'Color', 'w');
    
    subplot(3,2,3);
    tit = sprintf('$t = $ %.2f\n%d training data\n', t_star(i), N0);
    bar(x_star, u_star(:,i),'m');
    xlabel('$x$')
    ylabel('$u(t,x)$')
    axis tight
    set(gca,'FontSize',14);
    set(gcf, 'Color', 'w');
    title(tit);
    
    subplot(3,2,4);
    tit = sprintf('$t = $ %.2f\n%d training data\n', t_star(i+1), N1);
    bar(x_star, u_star(:,i+1),'m');
    xlabel('$x$')
    ylabel('$u(t,x)$')
    axis tight
    set(gca,'FontSize',14);
    set(gcf, 'Color', 'w');
    title(tit);
    
    subplot(3,2,5:6);

    s = '$\begin{tabular}{|c|c|}';
    s = strcat(s, ' \hline');
    s = strcat(s, ' Correct PDE & $u_t - 0.5 u_{xx} = 0$ \\');
    s = strcat(s, ' \hline');
    s = strcat(s, ' Identified PDE & $u_t - ', sprintf('%.3f',params(1)), '\mathcal{D}_{-\infty, x}^{', sprintf('%.3f',params(2)), '}u = 0$ \\');
    s = strcat(s, ' \hline');
    s = strcat(s, ' \end{tabular}$');    
    text(0.25,0.8,s,'interpreter','latex','FontSize',18)
    axis off
    
    if save_plt == 1
        export_fig ../Figures/Fractional.png -r300
    end
    
    drawnow();
end

end

function [t_star, x_star, u_star, pos] = gen_data(N,dt,m,n)

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
    
    f = figure();
    clf
    hold
    for i = 1:m
        h = histogram(P(:,i), bins, 'Normalization', 'pdf');
        u_star(:,i) = h.Values;
    end
    
    close(f)

end
