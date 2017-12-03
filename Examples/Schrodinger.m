% @author: Maziar Raissi

function Schrodinger()

clc; close all;

plt = 1;
plt_pred = 0;
save_plt = 1;

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

set(0,'defaulttextinterpreter','latex')

%% Load Data
load('../Data/nls.mat', 'usol', 't', 'x')
u_star = real(usol)'; % 501x512
v_star = imag(usol)'; % 501x512
t_star = t; % 501x1
x_star = x';   % 512x1
N_star = size(x_star,1);
nsteps = size(t_star,1)-1;
    
%% Setup
N0 = 49;
N1 = 51;
%% Clean Data
rng('default')
i = randi(nsteps);
dt = t_star(i+1) - t_star(i);
    
idx0 = randsample(N_star, N0);
x0 = x_star(idx0,:);
U0 = [u_star(idx0,i) v_star(idx0,i)];
    
idx1 = randsample(N_star,N1);
x1 = x_star(idx1,:);
U1 = [u_star(idx1,i+1) v_star(idx1,i+1)];

hyp = [log([1.0 1.0]) log([1.0 1.0]) 0.0 0.0 -4.0];
model = HPM(x1, U1, x0, U0, dt, hyp);
model = model.train(5000);
    
hyp = model.hyp;
params = hyp(5:6);
    
[pred_n_star, var_n_star] = model.predict(x_star);
pred_U_star = reshape(pred_n_star,N_star,2);
var_n_star = abs(diag(var_n_star));
var_U_star = reshape(var_n_star,N_star,2);
    
error_u = norm(pred_U_star(:,1) - u_star(:,i+1))/norm(u_star(:,i+1));
error_v = norm(pred_U_star(:,2) - v_star(:,i+1))/norm(v_star(:,i+1));
    
fprintf(1,'=========================\n');
fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error_u = %.2e, Error_v = %.2e\n\n', i, ...
    t_star(i+1), model.NLML, error_u, error_v);
       
str = sprintf('%.4f  ', params);
fprintf('Parameters: %s\n\n', str)
fprintf(1,'=========================\n\n');

if plt_pred == 1
    figure();
    subplot(2,1,1)
    plot_prediction_1D(x_star, u_star(:,i+1), pred_U_star(:,1), var_U_star(:,1), ...
        '$x$', '$u(t,x)$', 'Prediction (clean data)');
    
    subplot(2,1,2)
    plot_prediction_1D(x_star, v_star(:,i+1), pred_U_star(:,2), var_U_star(:,2), ...
        '$x$', '$v(t,x)$', 'Prediction (clean data)');
    
    drawnow;
end

%% Noisy Data
noise = 0.01;
U0(:,1) = U0(:,1) + noise*std(U0(:,1))*randn(size(U0(:,1)));
U0(:,2) = U0(:,2) + noise*std(U0(:,2))*randn(size(U0(:,2)));
U1(:,1) = U1(:,1) + noise*std(U1(:,1))*randn(size(U1(:,1)));
U1(:,2) = U1(:,2) + noise*std(U1(:,2))*randn(size(U1(:,2)));

hyp = [log([1.0 1.0]) log([1.0 1.0]) 0.0 0.0 -4.0];
model = HPM(x1, U1, x0, U0, dt, hyp);
model = model.train(5000);
    
hyp = model.hyp;
params_noise = hyp(5:6);
    
[pred_n_star, var_n_star] = model.predict(x_star);
pred_U_star = reshape(pred_n_star,N_star,2);
var_n_star = abs(diag(var_n_star));
var_U_star = reshape(var_n_star,N_star,2);
    
error_u = norm(pred_U_star(:,1) - u_star(:,i+1))/norm(u_star(:,i+1));
error_v = norm(pred_U_star(:,2) - v_star(:,i+1))/norm(v_star(:,i+1));
    
fprintf(1,'=========================\n');
fprintf(1,'Step: %d, Time = %.2f\n\nNLML = %.2f, Error_u = %.2e, Error_v = %.2e\n\n', i, ...
    t_star(i+1), model.NLML, error_u, error_v);
       
str = sprintf('%.4f  ', params_noise);
fprintf('Parameters: %s\n\n', str)
fprintf(1,'=========================\n\n');

if plt_pred == 1
    figure();
    subplot(2,1,1)
    plot_prediction_1D(x_star, u_star(:,i+1), pred_U_star(:,1), var_U_star(:,1), ...
        '$x$', '$u(t,x)$', 'Prediction (clean data)');
    
    subplot(2,1,2)
    plot_prediction_1D(x_star, v_star(:,i+1), pred_U_star(:,2), var_U_star(:,2), ...
        '$x$', '$v(t,x)$', 'Prediction (clean data)');
    
    drawnow;
end

%% Plot Results

if plt == 1
    fig = figure();
    set(fig,'units','normalized','outerposition',[0 0 1 0.75])
    subplot(5,2,1:2)
    plot_surface(t_star, x_star, u_star, '', '$x$', '$u(t,x)$');
    
    hold on
    plot3([t_star(i) t_star(i)],get(gca,'ylim'),[10 10],'k','LineWidth',2)
    plot3([t_star(i+1) t_star(i+1)],get(gca,'ylim'),[10 10],'k','LineWidth',2)
    
    subplot(5,2,3:4)
    plot_surface(t_star, x_star, v_star, '$t$', '$x$', '$v(t,x)$');
    
    hold on
    plot3([t_star(i) t_star(i)],get(gca,'ylim'),[10 10],'k','LineWidth',2)
    plot3([t_star(i+1) t_star(i+1)],get(gca,'ylim'),[10 10],'k','LineWidth',2)
    
    subplot(5,2,5);
    tit = sprintf('$t = $ %.5f\n%d training data\n', t_star(i), N0);
    plot_data_1D(x_star, u_star(:,i), x0, U0(:,1), '$x$', '$u(t,x)$', tit);
    
    subplot(5,2,7);
    plot_data_1D(x_star, v_star(:,i), x0, U0(:,2), '$x$', '$v(t,x)$', '');
    
    subplot(5,2,6);
    tit = sprintf('$t = $ %.5f\n%d training data\n', t_star(i+1), N1);
    plot_data_1D(x_star, u_star(:,i+1), x1, U1(:,1), '$x$', '$u(t,x)$', tit);
    
    subplot(5,2,8);
    plot_data_1D(x_star, v_star(:,i+1), x1, U1(:,2), '$x$', '$v(t,x)$', '');
    
    subplot(5,2,9:10);
    s1 = '$\begin{tabular}{ |c|c| }  \hline Correct PDE & $i h_t + 0.5 h_{xx} + |h|^2 h = 0$ \\  \hline Identified PDE (clean data) &';
    s2 = sprintf('$i h_t + %.3f h_{xx} + %.3f |h|^2 h = 0$', params(1), params(2));
    s3 = ' \\  \hline Identified PDE (1\% noise) &';
    s4 = sprintf('$i h_t + %.3f h_{xx} + %.3f |h|^2 h = 0$', params_noise(1), params_noise(2));
    s5 = ' \\  \hline \end{tabular}$';
    s = strcat(s1,s2,s3,s4,s5);
    text(0.1,0.8,s,'interpreter','latex','FontSize',18)
    axis off
    
    if save_plt == 1
        export_fig ../Figures/Schrodinger.png -r300
    end
    
    drawnow();
end

end