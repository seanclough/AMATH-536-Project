function plot_means(sheet)
    colors = ["#000000","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
    av_matrix = zeros(21,8);
    for i=1:8
        for j=1:21
            av_matrix(j,i) = mean(sheet(:,(7+21*(i-1)+j)));
        end
    end
    figure
    hold on
    for i=1:8
        scatter(0:200:4000,av_matrix(:,i),"MarkerEdgeColor",colors(i))
    end
end
