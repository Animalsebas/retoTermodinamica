% Constantes para el nitrógeno
a_N2 = 1.408; 
b_N2 = 0.039; 

n = 1; 
R = 8.314; 

V = linspace(0.001, 1, 10000); 

figure;

T = [0.5, 0.7, 0.9,  1.1, 1.3]; 
colors = hsv(length(T)); 

%for i = 1:length(T)
 %   P2 = (n * R * T(i)) ./ (V);
  %  plot(V, P2, '-');
   % hold on;
%end
hold on
for i = 1:length(T)
    P = (n * R * T(i)) ./ (V - n * b_N2) - a_N2 * n^2 ./ V.^2;
    plot(V, P, '.');
    hold on;
end

xlabel('Volumen');
ylabel('Presión');
title('Curvas Isotermas para el Nitrógeno');
legend('T = 0.5 K','T = 0.7 k','T = 0.9 K','T = 1.1 K', 'T = 1.3 K');


%hLegend = findobj(gcf, 'Type', 'Legend');
%set(hLegend, 'FontSize', 6); 

grid on;

% Marcador de punto de inflexión
Vc = 0.089; % Volumen molar crítico para el nitrógeno
Pc = 3.4e6; % Presión crítica para el nitrógeno
plot(0.117, 0.041, 'ro',  'DisplayName', 'Punto de Inflexión (Tc)');
ylim([-200, 400]);
xlim([0, 0.3]);