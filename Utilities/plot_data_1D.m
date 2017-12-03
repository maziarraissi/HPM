% @author: Maziar Raissi

function plot_data_1D(x_star, f_star, x_f, y_f, xlab, ylab, tit)

hold
plot(x_star, f_star, 'b', 'LineWidth', 2);
plot(x_f, y_f, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(xlab)
ylabel(ylab)
axis tight
set(gca,'FontSize',14);
set(gcf, 'Color', 'w');

title(tit);