function PlotMatrix(t,m)
    if size(m,1) == size (t)
        for i=1:size(m,2)
            plot(t,m(:,i));
            hold on;
        end
    else
        for i=1:size(m,1)
            plot(t,m(i,:));
            hold on;
        end
    end
    grid on;
end