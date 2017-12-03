% @author: Maziar Raissi

function Navier_Stokes()

clc; close all;

plt = 1;
plt_pred = 0;
save_plt = 1;

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

set(0,'defaulttextinterpreter','latex')

%% Load Data
load('../Data/cylinder_fine.mat', 't_star', 'X_star', 'U_star', 'w_star')
% load('../Data/cylinder_finer.mat', 't_star', 'X_star', 'U_star', 'w_star')
% load('../Data/cylinder_finest.mat', 't_star', 'X_star', 'U_star', 'w_star')
N_star = size(X_star,1);
nsteps = size(t_star,1) - 1;
    
%% Setup
N0 = 251;
N1 = 249;
%% Clean Data
rng('default')
i = 10;
% i = 19;
% i = 37;
% i = randi(nsteps);
dt = t_star(i+1) - t_star(i);
    
idx0 = randsample(N_star, N0);
X0 = X_star(idx0,:);
U0 = U_star(idx0,:,i);
    
idx1 = randsample(N_star,N1);
X1 = X_star(idx1,:);
U1 = U_star(idx1,:,i+1);

hyp = [log([1.0 1.0 1.0]) 0.0 0.0 log([1.0 1.0 1.0]) 0.0];
model = HPM(X1, U1, X0, U0, dt, hyp);
model = model.train(5000);
    
hyp = model.hyp;
params = [hyp(4) exp(hyp(5))];
    
[pred_n_star, var_n_star] = model.predict(X_star);
pred_U_star = reshape(pred_n_star,N_star,2);
var_n_star = abs(diag(var_n_star));
    
error_u = norm(pred_U_star(:,1) - U_star(:,1,i+1))/norm(U_star(:,1,i+1));
error_v = norm(pred_U_star(:,2) - U_star(:,2,i+1))/norm(U_star(:,2,i+1));
    
fprintf(1,'=========================\n');
fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error_u = %.2e, Error_v = %.2e\n\n', i, ...
    t_star(i+1), model.NLML, error_u, error_v);

str = sprintf('%.4f  ', params);
fprintf('Parameters: %s\n\n', str)
fprintf(1,'=========================\n\n');
    
if plt_pred == 1

    fig = figure();
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf

    subplot(2,2,1)
    plot_surface_griddata(X_star, U_star(:,1,i+1))

    subplot(2,2,2);
    plot_surface_griddata(X_star, U_star(:,2,i+1))

    subplot(2,2,3)
    plot_surface_griddata(X_star, pred_U_star(:,1))

    subplot(2,2,4);
    plot_surface_griddata(X_star, pred_U_star(:,2))       

    drawnow()
end

%% Noisy Data
noise = 0.01;

U0(:,1) = U0(:,1) + noise*std(U0(:,1))*randn(size(U0(:,1)));
U0(:,2) = U0(:,2) + noise*std(U0(:,2))*randn(size(U0(:,2)));
U1(:,1) = U1(:,1) + noise*std(U1(:,1))*randn(size(U1(:,1)));
U1(:,2) = U1(:,2) + noise*std(U1(:,2))*randn(size(U1(:,2)));

hyp = [log([1.0 1.0 1.0]) 0.0 0.0 log([1.0 1.0 1.0]) 0.0];
model = HPM(X1, U1, X0, U0, dt, hyp);
model = model.train(5000);
    
hyp = model.hyp;
params_noise = [hyp(4) exp(hyp(5))];
    
[pred_n_star, var_n_star] = model.predict(X_star);
pred_U_star = reshape(pred_n_star,N_star,2);
var_n_star = abs(diag(var_n_star));
    
error_u = norm(pred_U_star(:,1) - U_star(:,1,i+1))/norm(U_star(:,1,i+1));
error_v = norm(pred_U_star(:,2) - U_star(:,2,i+1))/norm(U_star(:,2,i+1));
    
fprintf(1,'=========================\n');
fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error_u = %.2e, Error_v = %.2e\n\n', i, ...
    t_star(i+1), model.NLML, error_u, error_v);

str = sprintf('%.4f  ', params_noise);
fprintf('Parameters: %s\n\n', str)
fprintf(1,'=========================\n\n');
    
if plt_pred == 1

    fig = figure();
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf

    subplot(2,2,1)
    plot_surface_griddata(X_star, U_star(:,1,i+1))

    subplot(2,2,2);
    plot_surface_griddata(X_star, U_star(:,2,i+1))

    subplot(2,2,3)
    plot_surface_griddata(X_star, pred_U_star(:,1))

    subplot(2,2,4);
    plot_surface_griddata(X_star, pred_U_star(:,2))       

    drawnow()
end

%% Plot Results

load('../Data/cylinder_vorticity.mat', 'XX', 'YY', 'WW')
if plt == 1
    
    fig = figure();
    set(fig,'units','normalized','outerposition',[0 0 1 .75])
    subplot(5,2,1:4)    
    tit = sprintf('$w(t,x,y)$, $t = $ %.2f', t_star(1));
    plot_surface(XX(1,:)', YY(:,1), WW, '$x$', '$y$', tit)

    xmin = 100;
    xmax = 425;
    ymin = 15;
    ymax = 185;

    hold on
    plot3([XX(1,xmin) XX(1,xmax)], [YY(ymin,1) YY(ymin,1)] ,[10 10],'k','LineWidth',3)
    plot3([XX(1,xmin) XX(1,xmax)], [YY(ymax,1) YY(ymax,1)] ,[10 10],'k','LineWidth',3)
    plot3([XX(1,xmin) XX(1,xmin)], [YY(ymin,1) YY(ymax,1)] ,[10 10],'k','LineWidth',3)
    plot3([XX(1,xmax) XX(1,xmax)], [YY(ymin,1) YY(ymax,1)] ,[10 10],'k','LineWidth',3)
    
    subplot(5,2,5)
    tit = sprintf('$u(t,x,y)$, $t = $ %.2f, %d training data', t_star(i), N0);
    plot_surface_griddata(X_star, U_star(:,1,i), '$x$', '$y$', tit)
    hold on
    plot3(X0(:,1), X0(:,2), 10*ones(size(X0(:,1))), 'k+', 'LineWidth',2)
    
    subplot(5,2,6)
    tit = sprintf('$u(t,x,y)$, $t = $ %.2f, %d training data', t_star(i+1), N1);
    plot_surface_griddata(X_star, U_star(:,1,i+1), '$x$', '$y$', tit)
    hold on
    plot3(X1(:,1), X1(:,2), 10*ones(size(X1(:,1))), 'k+', 'LineWidth',2)
    
    subplot(5,2,7)
    tit = sprintf('$v(t,x,y)$', t_star(i), N0);
    plot_surface_griddata(X_star, U_star(:,2,i), '$x$', '$y$', tit)
    hold on
    plot3(X0(:,1), X0(:,2), 10*ones(size(X0(:,1))), 'k+', 'LineWidth',2)
    
    subplot(5,2,8)
    tit = sprintf('$v(t,x,y)$', t_star(i+1), N1);
    plot_surface_griddata(X_star, U_star(:,2,i+1), '$x$', '$y$', tit)
    hold on
    plot3(X1(:,1), X1(:,2), 10*ones(size(X1(:,1))), 'k+', 'LineWidth',2)
    
    subplot(5,2,9:10);
    s = '$\begin{tabular}{|c|c|}';
    s = strcat(s, ' \hline');
    s = strcat(s, ' Correct PDE & $\begin{array}{c}');
    s = strcat(s, ' u_t + (u u_x + v u_y) = -p_x + 0.01 (u_{xx} + u_{yy})\\');
    s = strcat(s, ' v_t + (u v_x + v v_y) = -p_y + 0.01 (v_{xx} + v_{yy})');
    s = strcat(s, ' \end{array}$ \\ ');
    s = strcat(s, ' \hline');
    s = strcat(s, ' Identified PDE (clean data) & $\begin{array}{c}');
    s = strcat(s, sprintf(' u_t + %.3f (u u_x + v u_y) = -p_x + %.5f (u_{xx} + u_{yy})', params(1), params(2)));
    s = strcat(s, ' \\');
    s = strcat(s, sprintf(' v_t + %.3f (u v_x + v v_y) = -p_y + %.5f (v_{xx} + v_{yy})', params(1), params(2)));
    s = strcat(s, ' \end{array}$ \\ ');
    s = strcat(s, ' \hline');
    s = strcat(s, ' Identified PDE (1\% noise) & $\begin{array}{c}');
    s = strcat(s, sprintf(' u_t + %.3f (u u_x + v u_y) = -p_x + %.5f (u_{xx} + u_{yy})', params_noise(1), params_noise(2)));
    s = strcat(s, ' \\');
    s = strcat(s, sprintf(' v_t + %.3f (u v_x + v v_y) = -p_y + %.5f (v_{xx} + v_{yy})', params_noise(1), params_noise(2)));
    s = strcat(s, ' \end{array}$ \\ ');     
    s = strcat(s, ' \hline');
    s = strcat(s, ' \end{tabular}$');
    text(0.0,0.5,s,'interpreter','latex','FontSize',18)
    axis off
    
    if save_plt == 1
        export_fig ../Figures/NavierStokes.png -r300
    end
    
    drawnow();
end

end