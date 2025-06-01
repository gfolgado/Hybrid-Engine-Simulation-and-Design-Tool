function Plotting(data, Pc)
%PLOTTING From CEA data, presents plots for user analysis
%   Allows user to select which plots to visualize

g = 9.8067;     % Acceleration due to gravity [m/s^2]

% Prompt the user for which graphs they wish to visualize
fprintf('Select which graphs to visualize:\n');
fprintf('1. OF vs Optimum Expansion Ratio\n');
fprintf('2. OF vs Characteristic Velocity\n');
fprintf('3. OF vs Thrust Coefficient\n');
fprintf('4. OF vs Vacuum Isp\n');
fprintf('5. OF vs Isp\n');
selected_graphs = input('Enter the numbers of the graphs to visualize (e.g., [1 2 3]): ');

% Initialize the color palette for plotting
color_palette = parula(size(data, 3));

% Iterate over selected graphs and plot each one
for graph_index = selected_graphs
    switch graph_index
        case 1
            % OF vs epsillon
            figure;
            hold on;
            for i = 1:size(data, 3)
                plot(data(:, 1, i), data(:, 2, i), '-o', 'Color', color_palette(i, :));
            end
            hold off;
            xlabel('OF ratio');
            ylabel('Expansion ratio');
            legend(cellstr(strcat(num2str(Pc'), ' [bar]')), 'Location', 'best');
            grid on;
            
        case 2
            % OF vs c*
            figure;
            hold on;
            for i = 1:size(data, 3)
                plot(data(:, 1, i), data(:, 3, i), '-o', 'Color', color_palette(i, :));
            end
            hold off;
            xlabel('OF ratio');
            ylabel('c*');
            legend(cellstr(strcat(num2str(Pc'), ' [bar]')), 'Location', 'best');
            grid on;
            
        case 3
            % OF vs Cf
            figure;
            hold on;
            for i = 1:size(data, 3)
                plot(data(:, 1, i), data(:, 4, i), '-o', 'Color', color_palette(i, :));
            end
            hold off;
            xlabel('OF ratio');
            ylabel('Cf');
            legend(cellstr(strcat(num2str(Pc'), ' [bar]')), 'Location', 'best');
            grid on;
            
        case 4
            % OF vs Vacuum Isp
            figure;
            hold on;
            for i = 1:size(data, 3)
                plot(data(:, 1, i), data(:, 5, i)/g, '-o', 'Color', color_palette(i, :));
            end
            hold off;
            xlabel('OF ratio');
            ylabel('Vacuum Isp');
            legend(cellstr(strcat(num2str(Pc'), ' [bar]')), 'Location', 'best');
            grid on;            
            
        case 5
            % OF vs Isp
            figure;
            hold on;
            for i = 1:size(data, 3)
                plot(data(:, 1, i), data(:, 6, i)/g, '-o', 'Color', color_palette(i, :));
            end
            hold off;
            xlabel('OF ratio');
            ylabel('Isp');
            legend(cellstr(strcat(num2str(Pc'), ' [bar]')), 'Location', 'best');
            grid on;            
        otherwise
            fprintf('Invalid selection. Skipping graph %d.\n', graph_index);
    end
end
end

