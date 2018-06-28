function DrawPath(idx, data, clr, lw)
hold on
    plot3(data(idx, 1), data(idx, 2),data(idx, 3), clr, 'LineWidth', lw);
hold off
