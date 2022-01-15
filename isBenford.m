%% This function extracts leading digits from a dataset, outputs the 
%distribution of each digit in percentage of the whole dataset,
%and plots it
function [benford, distrib] = isBenford(X, varargin)

hlp_varargin2struct(varargin, ... 
    {'plotting','plot'},'on', ...
    {'printing','print'},'on');

% disp('processing...');

X = X(~isnan(X));               %Remove NaN values

for i = 1:length(X)
    string = sprintf('%0.5e', abs(X(i)));       %Convert values to string
    firstDigit(i,:) = str2double(string(1));    %Keep only leading digit
end
firstDigit = firstDigit(firstDigit ~= 0);       %Remove 0s

%Get % for each value
one = sum(firstDigit == 1) / length(firstDigit);                     
two = sum(firstDigit == 2) / length(firstDigit);                     
three = sum(firstDigit == 3) / length(firstDigit);                     
four = sum(firstDigit == 4) / length(firstDigit);                     
five = sum(firstDigit == 5) / length(firstDigit);                     
six = sum(firstDigit == 6) / length(firstDigit);                     
seven = sum(firstDigit == 7) / length(firstDigit);                     
eight = sum(firstDigit == 8) / length(firstDigit);                     
nine = sum(firstDigit == 9) / length(firstDigit);

distrib = [one, two, three, four, five, six, seven, eight, nine];

benford = log10(1 + 1 ./ [1 2 3 4 5 6 7 8 9]);

%Plot
if strcmp(plotting,'on')
    figure; bar(distrib); hold on; plot(benford, 'r', 'LineWidth',2); legend('data', 'Benford law')
    grid on;
    xlabel('Leading digit');
    ylabel('Percent');
end

%Calculate and display ratios
if strcmp(printing,'on')
    fprintf(['Benford ratios: ' num2str(round(benford(1)/benford(2),1)) '\t' num2str(round(benford(2)/benford(3),1)) '\t' num2str(round(benford(3)/benford(4),1)) '\t' num2str(round(benford(4)/benford(5),1)) '\t' num2str(round(benford(5)/benford(6),1)) '\t' num2str(round(benford(6)/benford(7),1)) '\t' num2str(round(benford(7)/benford(8),1)) '\t' num2str(round(benford(8)/benford(9),1)) '\n']);
    fprintf(['Data ratios:    ' num2str(round(one/two,1)) '\t' num2str(round(two/three,1)) '\t' num2str(round(three/four,1)) '\t' num2str(round(four/five,1)) '\t' num2str(round(five/six,1)) '\t' num2str(round(six/seven,1)) '\t' num2str(round(seven/eight,1)) '\t' num2str(round(eight/nine,1)) '\n'])
end

end