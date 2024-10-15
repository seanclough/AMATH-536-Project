%% plot means
%plot_means(complete_sheet);
%plot_theo()
plot_means(data_sheet);
plot_theo();
xlabel('time (days)');
ylabel('cells');
legend_titles = {'stem'};
for i=1:7
    legend_titles{end+1}=string(i);
end
legend(legend_titles,'Location','northwest');
legend('boxoff');

%% plot variances
av_matrix = zeros(21,8);
for i=1:8
    for j=1:21
        sum = 0;
        for k=1:length(sheet)
            %compartments start in column 8 and 
            %there are 21 columns per compartment
            sum = sum + sheet(k,(7+21*(i-1)+j));
        end
        %av_matrix(j,i) = sum/length(sheet);
        av_matrix(j,i) = mean(sheet(:,(7+21*(i-1)+j)));
    end
end
var_matrix = zeros(21,8);
for i=1:8
    for j=1:21
        sum = 0;
        for k=1:length(sheet)
            %(x - x_bar)^2
            sum = (sheet(k,(7+21*(i-1)+j))-av_matrix(j,i))^2;
        end
        %using unbiased estimator
        %var_matrix(j,i) = sum/(length(sheet)-1);
        var_matrix(j,i) = var(sheet(:,(7+21*(i-1)+j)));
    end
end
cv_matrix = zeros(21,8);
for i=1:8
    for j=1:21
        %CV = stdev/x_bar
        %matlab must be coded to skip dividing by 0...
        cv_matrix(j,i) = sqrt(var_matrix(j,i))/av_matrix(j,i);
    end
end
tiledlayout(2,4)
colors = ["#000000","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
for i=1:8
    nexttile
    scatter(0:200:4000,var_matrix(:,i),"MarkerEdgeColor",colors(i))
    %scatter(0:200:4000,cv_matrix(:,i),"MarkerEdgeColor",colors(i))
    title(append("Compartment ",string(i-1)));
end

%% testing area
count_extinctions = 0;
for i=1:length(sheet)
    if sheet(i,28)<1
        count_extinctions = count_extinctions + 1;
    end
end
count_extinctions


%{
parfor i=1:100000
    ssa_project
end 
%}

% times = 10:10:100;
% result = ssa_try(times,10,1,1,100000,false,false);
% result;
% av_result=zeros(1,length(result(:,1)));
% for i=1:length(av_result)
%    av_result(i)=sum(result(i,:))/length(result(i,:));
% end
% av_result
% %av_result = [92,308,764,1781,5569,14820,40111,108097,294187,799294]
% figure()
% scatter(times, av_result)
% %set(gca, "YScale", "log");
% xlabel("time (s)");
% ylabel("Average Number of Cells");

% close all
% figure();
% hold on;
% xlabel("time (s)");
% ylabel("cells");
% for i=1:1
%     [finalx,t,x]=onerun(100);
%     stairs(t,x)
%     finalx
% end
% %set(gca, "YScale", "log");
% set(gca, "XLim", [0,t(length(t))])

%x = simulation(100,10000);
%x