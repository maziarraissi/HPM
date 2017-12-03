% @author: Maziar Raissi

function plot_surface(x_star, y_star, Z_star, xlab, ylab, tit)

[X_star, Y_star] = meshgrid(x_star, y_star);
surface(X_star, Y_star, Z_star);
axis tight
colormap jet
shading interp
colorbar
xlabel(xlab);
ylabel(ylab);
title(tit);
set(gca,'FontSize',14);
set(gcf, 'Color', 'w');