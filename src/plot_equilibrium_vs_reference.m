function plot_equilibrium_vs_reference(equilibrium,reference,row,col,~)
% Optional visualization used only by verbose reference diagnostics.
N = row*col;
figure; hold on;
plot(equilibrium(1:N),equilibrium(N+1:2*N),'o','DisplayName','Rest configuration');
plot(reference(1:N),reference(N+1:2*N),'x','DisplayName','Loaded reference');
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); legend('show');
end
