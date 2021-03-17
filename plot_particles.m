function plot_particles(as,a)
% Plot the particles in state space when k=3
figure;
plot3(as(:,1),as(:,2),as(:,3),'o')
xlim([-.45,.45]);ylim([-.3,.3]);zlim([-.25,.25])
xlabel('a_1');ylabel('a_2');zlabel('a_3');
box on;hold on;
plot3(a(1),a(2),a(3),'r*','linewidth',10)
