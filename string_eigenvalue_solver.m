function [ naturalFrequencies, naturalModes ] = ...
string_eigenvalue_solver( masses, F, h )
% string_eigenvalue_solver Solves simple eigenvalue problem of a string
%   Solves simple eigenvalue problem of a string with user-given masses in
%   a column or row vector, constant tension force F on the string and
%   separation of masses h

% Get dimension of arrays by size of masses vector
% Getting the maximum of the dimensions will allow the mass vector
% to be a column vector or row vector
n = max(size(masses,1), size(masses,2));

% --- Generate A ---
% 2s on main diagonal
A = 2*eye(n);

% -1s on second diagonals
diagOnes = -1 * ones(n-1,1);
A = A + diag(diagOnes,1) + diag(diagOnes,-1);

% --- Generate D ---
% Diagonal elements
diags = F ./ (h * masses);
D = diag(diags);

% Get eigenvalues and eigenvectors (natrual frequencies) of D * A
[naturalModes, eigenVals] = eig(D * A);

% Natural frequencies are square roots of eigenvalues times -i
naturalFrequencies = -1i * sqrt(eigenVals);

% --- Plotting ---
% Range of numbers for plotting
maxVal = h * n;
x = h:h:maxVal;

% Create new figure
plotFig = figure;
plotAx = axes;
hold(plotAx);

% Markers for plotting
markers = {'--+',':o','-.*','--x',':s','-.d', ...
    '--^',':v','-.>','--<',':p','-.h'};

% Plot eigenmodes with frequencies
for i=1:size(naturalModes, 2)
    plot(plotAx, x, naturalModes(:,i), markers{mod(i,numel(markers))}, ...
        'DisplayName', sprintf('Mode %d, Frequency = %f', ...
        i, eigenVals(i,i)), 'Linewidth', 2);
end

% Set limits of plot
xlim(plotAx, [min(x) - 1, max(x) + 1]);

% Show legend
plotLegend = legend('-DynamicLegend');
set(plotLegend, 'FontSize', 30);
set(plotLegend, 'Location', 'southoutside');

% Toggle off hold on axes
hold(plotAx);
end

