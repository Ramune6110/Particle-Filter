clear;
close all;
clc;
% reading file
fileID = fopen('data.txt', 'r');
formatSpec = '%f';
size = [5 5000];
data = fscanf(fileID, formatSpec, size);
% graph
figure(1);
set(groot, 'DefaultTextInterpreter', 'Latex');
set(groot, 'DefaultLegendInterpreter', 'Latex');

figure(1);
%Plot Result
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(data(2, :), data(3, :),'-.b'); hold on;
plot(data(4, :), data(5, :),'r'); hold on;
title('Particle Filter', 'fontsize', 16, 'fontname', 'times');
xlabel('X (m)', 'fontsize', 16, 'fontname', 'times');
ylabel('Y (m)', 'fontsize', 16, 'fontname', 'times');
legend('Ground Truth','PF');
grid on;
axis equal;

% Auto save graph
saveas(gcf, 'Particle_Filter.png');