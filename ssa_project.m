function ssa_project()
% Simulate stem cell + first 6 non stem cell populations
import Gillespie.*
%close all
global n
%% Reaction network:
%   1. stem cell division:                    0 --N_0_div--> N_0
%   2. stem cell differentiation            N_0 --N_0_dif--> N_1
%   3. ith compartment division:              0 --N_i_div--> N_i
%   4. ith compartment differentiation:     N_i --N_i_dif--> N_{i+1}
%   5. ith compartment death:               N_i --N_i_dth--> 0

%% Initial state
tspan = [0, 4001]; %days
notetimes = zeros(1,21); %times to be recorded
for i=1:21
    notetimes(i) = (i-1)*200;
end 
n = 7; %number of non stem cell compartments
x0 = zeros(1,n+1); 
x0(1) = 400; %stem cell population
%x0(2) = 1; %one mutant cell in 1st non stem cell compartment

%% Rate constants
dif_rate = 0.85;
dth_rate = 0;
p.N_0_div = 0.5;
p.N_0_dif = 0.5;
p.N_i_div = 1 - dif_rate - dth_rate;
p.N_i_dif = dif_rate;
p.N_i_dth = dth_rate;

prolif_rate = 1.26;
p.base_rate = 1/365;
p.prolif_rate = zeros(1,n+1);
for i=1:n
    p.prolif_rate(i) = p.base_rate*(prolif_rate)^i;
end
%precalculation should hopefully speed up the code

%% Specify reaction network
pfun = @propensities;
stoich_matrix = zeros(3*n+2,n+1);
stoich_matrix(1,1)=1; stoich_matrix(2,1)=-1; stoich_matrix(2,2)=1;
for i=1:n
    stoich_matrix(3*i,i+1)=1;
    stoich_matrix(3*i+1,i+1)=-1; stoich_matrix(3*i+1,i+2)=2;
    stoich_matrix(3*i+2,i+1)=-1;
end
stoich_matrix=stoich_matrix(1:3*n+2,1:n+1); %assumes there is another compartment after last one

%% Run simulation
[t,x,r] = directMethod(stoich_matrix, pfun, tspan, x0, p, notetimes);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
%{
figure();
stairs(t,x); 
set(gca,'XLim',tspan);
set(gca, 'YScale', 'linear');
xlabel('time (days)');
ylabel('cells');
legend_titles = {'stem'};
for i=1:n
    legend_titles{end+1}=string(i);
end
legend(legend_titles,'Location','northwest');
legend('boxoff');
%}

%% Record Data
load("data_sheet.mat","data_sheet");
now_time = datetime("now");
data_row = zeros(1,175);
data_row(1) = 60*hour(now_time)+minute(now_time)+1/60*second(now_time);
data_row(2) = day(now_time,"dayofyear");
data_row(3) = year(now_time);
data_row(4) = p.N_i_dif;
data_row(5) = p.N_i_dth;
data_row(6) = prolif_rate;
data_row(7) = p.base_rate;
for i=1:(n+1)
    data_row(21*(i-1)+8:21*(i-1)+28) = r(i,:);
end
data_sheet = [data_sheet; data_row];
save("data_sheet.mat","data_sheet");
%% Testing zone
%{
test1 = [1,2,3,4,5,6,7];
test2 = propensities(test1,p);
test1
test2
%}

%{
test3 = zeros(1,n+1);
for i=1:n+1
    test3(i) = 400*(2*.85)^(i-2)/((1.26*0.7)^(i-1));
end
test3/1000
%}
end



function a = propensities(x, p)
% Return reaction propensities given current state x
    global n
    a = zeros(3*n+2,1);
    a(1) = p.N_0_div*p.base_rate*x(1);
    a(2) = p.N_0_dif*p.base_rate*x(1);
    for i=1:n
        a(3*i)   = p.N_i_div*p.prolif_rate(i)*x(i+1);
        a(3*i+1) = p.N_i_dif*p.prolif_rate(i)*x(i+1);
        a(3*1+2) = p.N_i_dth*p.prolif_rate(i)*x(i+1);
    end
end

