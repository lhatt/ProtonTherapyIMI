function plot_particles(as,a)
% Plot the particles in state space when k=3
%figure;
scatter3(as(:,1),as(:,2),as(:,3))
xlim([0,2]),ylim([0,2]),zlim([0,2])
xlabel('d_1');ylabel('d_2');zlabel('d_3');
box on;grid on;hold on;
scatter3(a(1),a(2),a(3),'r*','linewidth',10)
