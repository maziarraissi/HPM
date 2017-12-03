% @author: Maziar Raissi

function plot_prediction_1D(x_star, f_star, mean_star, var_star, ...
    xlab, ylab, tit)

color = [217,95,2]/255;

hold
plot(x_star, f_star, 'b', 'LineWidth', 3);
plot(x_star, mean_star,'r--','LineWidth',3);
[l,p] = boundedline(x_star, mean_star, 2.0*sqrt(var_star), ':', 'alpha','cmap', color);
outlinebounds(l,p);
xlabel(xlab)
ylabel(ylab)
axis tight
set(gca,'FontSize',14);
set(gcf, 'Color', 'w');

title(tit);