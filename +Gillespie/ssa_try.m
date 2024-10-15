function ssa_try()
% Simulate my homework
import Gillespie.*

%% Reaction network:
%   1. birth:       0       --kl--> cell
%   2. death:       cell    --km--> 0

%% Rate constants
p.kl = 1;%0.01;      
p.km = 1;%1;                     

%% Initial state
tspan = [0, 100]; %seconds
x0    = 10;     %cell inital pop

%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ 1    %birth
                 -1 ]; %pdeath

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
figure();
stairs(t,x); set(gca,'XLim',tspan);
xlabel('time (s)');
ylabel('cells');

end


function a = propensities_2state(x, p)
% Return reaction propensities given current state x
cell  = x(1);

a = [p.kl*cell;    %birth
     p.kd*cell];   %death
end
