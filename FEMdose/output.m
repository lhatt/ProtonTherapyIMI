%This script gives some examples for the output of interesting things for
%the report writeup.

%Plot the final mesh in the algorithm
% figure(1)
% triplot(mesh.elem_vertices, ...
%     mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
% axis equal
% title('final mesh');


%Plot the exact solution against the finite element approximation
figure(2)
% coord = [mesh.vertex_coordinates(:,1),mesh.vertex_coordinates(:,2)];
% for i = 1:length(mesh.vertex_coordinates(:,1))
%     ex(i) = feval(u_exact,coord(i,:));
% end
% subplot(1,2,1)
% trisurf(mesh.elem_vertices, ...
%     mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
%     ex);
% title('solution')
% subplot(1,2,2)
trisurf(mesh.elem_vertices, ...
    mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
    uh);
title('FE solution U')
xlabel('x')
ylabel('y')

%Plot the convergence rate of the H1 error and the estimator
% figure(3);
% xlim([0 10^5]);
% loglog(ndof,err_H1,'r -o','linewidth',2);
% hold on;
% xlabel('dim V_h');
% ylabel('value');
% loglog(ndof, est,'b -o','linewidth',2);
% str = ['EOC = ' num2str(EOC_H1(end))];
% text(ndof(end),err_H1(end),str,...
%     'VerticalAlignment','top','HorizontalAlignment','right')
% str = ['EOC = ' num2str(EOC_est(end))];
% text(ndof(end),est(end),str,...
%     'VerticalAlignment','bottom','HorizontalAlignment','left')
% legend('H^1 error','estimator');