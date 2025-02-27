% Prompt user to input the directory where data files are located
folder = uigetdir('', 'Select the folder containing data files');

% Prompt user for the sample name
prompt = {'Enter the sample name:'};
dlgtitle = 'Sample Name Input';
dims = [1 35];
definput = {'Sample1'};
sampleName = inputdlg(prompt, dlgtitle, dims, definput);

% Check if the user canceled the input dialog
if isempty(sampleName)
    disp('Sample name input cancelled. Script terminated.');
    return;
end

% Use the sample name for the bar plot title
sampleName = sampleName{1};  % Extract the input from the cell array
% Check if the folder selection was cancelled
if folder == 0
    disp('Folder selection cancelled. Script terminated.');
    return;
end

% List all files in the specified folder with the '.dat' extension
files = dir(fullfile(folder, '*.dat'));

% Filter files that contain 'trace' in their names (case insensitive)
traceFiles = files(contains(lower({files.name}), 'tr'));

% Check if no trace files were found
if isempty(traceFiles)
    disp('No files found with "trace" in their names. Script terminated.');
    return;
end

% Initialize arrays to store results
averages1 = zeros(length(traceFiles), 1);
averages2 = zeros(length(traceFiles), 1);
averages3 = zeros(length(traceFiles), 1);
fileNames = cell(length(traceFiles), 1);

% Loop through each trace file
for i = 1:length(traceFiles)
    % Construct the full file path
    filePath = fullfile(folder, traceFiles(i).name);
    
   
        % Read data from the file
        data = importdata(filePath);
        
       % Display the first five rows of the data
    if isstruct(data)  % Check if importdata returns a struct
        disp('First five rows of the data:');
        disp(data.data(1:min(5, size(data.data, 1)), :));
        secondColumn = data.data(:, 2);
    else
        disp('First five rows of the data:');
        disp(data(1:min(5, size(data, 1)), :));
        secondColumn = data(:, 2);
    end
        
        % Calculate the average of the first 10 values in the second column, excluding the first data point
        averages1(i) = mean(secondColumn(2:15));
        
        % Calculate the average from 50 to 100 values in the second column
        averages2(i) = mean(secondColumn(50:100));
        
        % Calculate the average from 450 to 500 values in the second column
        averages3(i) = mean(secondColumn(450:500));
        
        % Store the filename
        fileNames{i} = traceFiles(i).name;
        
        % Display results for each file (optional)
        fprintf('File: %s\n', traceFiles(i).name);
        fprintf('  Average of first 10 values in second column (excluding first): %.2f\n', averages1(i));
        fprintf('  Average from 50 to 100 in second column: %.2f\n', averages2(i));
        fprintf('  Average from 450 to 500 in second column: %.2f\n', averages3(i));
    %catch
        % Display error message if unable to read data
        fprintf('Error reading data from file: %s\n', traceFiles(i).name);
    end
%end

% Remove zero entries if any due to errors (assuming average values can't be exactly zero)
averages1 = averages1(averages1 ~= 0);
averages2 = averages2(averages2 ~= 0);
averages3 = averages3(averages3 ~= 0);

% Filter out values less than zero in averages3
averages3 = averages3(averages3 >= 0);

% Check if there are any valid averages to plot
if isempty(averages1) && isempty(averages2) && isempty(averages3)
    disp('No valid data found in files. Plot cannot be generated.');
    return;
end

% Calculate the mean value of the filtered averages3
meanAverages3 = mean(averages3);
fprintf('Mean of filtered averages3: %.2f\n', meanAverages3);

% Subtract the mean value of averages3 from averages1 and averages2
differences1 = averages1 - meanAverages3;
differences2 = averages2 - meanAverages3;

%filtering values less than 100 and larger than 700
differences1 = differences1(differences1 >= 100 & differences1 <= 700);
differences1 = differences1/300

% Specify bin size (bin width) for the histogram
binSize = 0.065;  % Adjust this value as desired

% Specify the x and y limits
xLimits = [0, 3];  % Adjust these values as needed
yLimits = [0, 2];    % Adjust these values as needed

% Plot histogram of differences1 with specified bin size
figure;
histogram(differences1, 'BinWidth', binSize, 'Normalization', 'pdf');
xlabel('Normalized intensity');
ylabel('Probability Density');
title(['Population Distribution of ', sampleName], 'FontSize', 14, 'FontWeight', 'bold');
xlim(xLimits);  % Set the x-axis limits
ylim(yLimits);  % Set the y-axis limits

% Add total number of difference values to the plot
numDifferences1 = length(differences1);
annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', ['n = ' num2str(numDifferences1)], 'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');

% Fit a Gaussian Mixture Model (GMM) to the differences1 data
%options = statset('MaxIter', 500);
%gm = fitgmdist(differences1, 2, 'Start', struct('mu', [220; 440], 'Sigma', cat(3, var(differences1), var(differences1)), 'ComponentProportion', [0.5 0.5]), 'Options', options);

% Define the options for the GMM fitting
options = statset('MaxIter', 1000, 'Display', 'final');

% Define initial parameters for the GMM
initialMu = [1; 2]; % Initial means
initialSigma = cat(3, var(differences1), var(differences1)); % Initial variances
initialProportions = [0.5, 0.5]; % Initial proportions

% Structure for initial parameters
initialParams = struct('mu', initialMu, 'Sigma', initialSigma, 'ComponentProportion', initialProportions);

% Fit the GMM with specified initial parameters
gm = fitgmdist(differences1, 2, 'Start', initialParams, 'Options', options);


% Get the x values for the fitted curve
x_values = linspace(min(differences1), max(differences1), 1000);

% Plot each component of the GMM separately
hold on;
component1 = gm.ComponentProportion(1) * normpdf(x_values, gm.mu(1), sqrt(gm.Sigma(1, 1, 1)));
component2 = gm.ComponentProportion(2) * normpdf(x_values, gm.mu(2), sqrt(gm.Sigma(1, 1, 2)));
plot(x_values, component1, 'b--', 'LineWidth', 1.5);
plot(x_values, component2, 'g--', 'LineWidth', 1.5);

% Add legend
%legend('Data', 'Component 1', 'Component 2');

% Optionally, save the plot as an image
outputPlot_diff1 = ['differences1_histogram_bin_' num2str(binSize) '.png'];  % Specify the output plot file name
saveas(gcf, outputPlot_diff1);
fprintf('Histogram plot saved as %s\n', outputPlot_diff1);

% Plot histogram of differences2 with specified bin size
%figure;
%histogram(differences2, 'BinWidth', binSize, 'Normalization', 'pdf');
%xlabel('Difference Value');
%ylabel('Probability Density');
%title('Histogram of Differences (Averages2 - Mean of Averages3)');

% Add total number of difference values to the plot
%numDifferences2 = length(differences2);
%annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', ['n = ' num2str(numDifferences2)], 'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');

% Optionally, save the plot as an image
%outputPlot_diff2 = ['differences2_histogram_bin_' num2str(binSize) '.png'];  % Specify the output plot file name
%saveas(gcf, outputPlot_diff2);
%fprintf('Histogram plot saved as %s\n', outputPlot_diff2);

% Calculate the area under each peak
area1 = gm.ComponentProportion(1);
area2 = gm.ComponentProportion(2);

%fprintf('Area under the first peak: %.2f\n', area1);
%fprintf('Area under the second peak: %.2f\n', area2);

% Specify the x and y limits
%yLimits = [0, 100];  % Adjust these values as needed

% Create a bar plot for area1 and area2
%figure;
%bar([area1*100, area2*100]);
%set(gca, 'XTickLabel', {'Monomer', 'Dimer'});
%ylabel('Percentage');
%title('Area under Each Peak');
%ylim(yLimits);  % Set the x-axis limits
%outputPlot_area = 'area_under_peaks.png';
%saveas(gcf, outputPlot_area);
%fprintf('Bar plot saved as %s\n', outputPlot_area);

%legend('Monomer = ' area1, 'Dimer = ' area2, 'n =' numDifferences2);
%annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', ['Monomer = ' num2str(area1*100),'% ', 'Dimer = ' num2str(area2*100),'% \n', 'n = ' num2str(numDifferences2),], 'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');
%annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', sprintf('Monomer = %.2f%%\nDimer = %.2f%%\nn = %d', area1*100, area2*100, numDifferences1), 'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');

% to find the standard deviation and error in the distribution
% Fit the GMM with specified initial parameters (as done previously)
% [Assuming gm has already been fitted with the specified parameters]

% Calculate the standard deviation for each component
stdDev1 = sqrt(gm.Sigma(1, 1, 1));
stdDev2 = sqrt(gm.Sigma(1, 1, 2));

% Calculate the standard error for each component
n1 = round(gm.ComponentProportion(1) * numDifferences1);
n2 = round(gm.ComponentProportion(2) * numDifferences1);
stdErr1 = stdDev1 / sqrt(n1);
stdErr2 = stdDev2 / sqrt(n2);

% Display the areas, standard deviations, and standard errors
fprintf('Area under first peak: %.4f, Std Dev: %.4f, Std Err: %.4f\n', area1, stdDev1, stdErr1);
fprintf('Area under second peak: %.4f, Std Dev: %.4f, Std Err: %.4f\n', area2, stdDev2, stdErr2);

% Specify the x and y limits
yLimits = [0, 100];  % Adjust these values as needed

% Create a bar plot to visualize the areas with error bars
figure;
bar([1, 2], [area1*100, area2*100]);
ylim(yLimits);  % Set the y-axis limits
hold on;
errorbar([1, 2], [area1*100, area2*100], [stdErr1*100, stdErr2*100], 'k.', 'LineWidth', 0.5);

% Customize the plot
set(gca, 'xticklabel', {'Monomer', 'Dimer'});
ylabel('Percentage');
title(['Percent oligomer of ', sampleName], 'FontSize', 14, 'FontWeight', 'bold');


annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', sprintf('Monomer = %.2f%%\nDimer = %.2f%%\nn = %d', area1*100, area2*100, numDifferences1), 'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'white', 'HorizontalAlignment', 'center');


% Optionally, save the plot as an image
outputPlot_with_error = ['differences1_histogram_bin_' num2str(binSize) '_with_gmm_error.png'];
saveas(gcf, outputPlot_with_error);
fprintf('Bar plot with standard error saved as %s\n', outputPlot_with_error);

% Save `differences1` values to a CSV file
outputFile_diff1 = fullfile(folder, ['differences1_', sampleName, '.csv']);
csvwrite(outputFile_diff1, differences1);
fprintf('Differences1 values saved to %s\n', outputFile_diff1);
