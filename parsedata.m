function QF = parsedata(data)
% Sarah Sundius, February 7, 2023

% Function to parse data files from Ellinor and extract "quantities final"
% Inputs:   data = data table
% Outputs:  QF = bacterial density

% Convert data table to a matrix
data = data{:,:};

% Get size
[r,c] = size(data);
numreps = data(end-1,1);
if isnan(numreps)
    numreps = data(end-2,1);
end

% Get initial condition
qf0 = data(end,4);

% Set up new matrix
QF = zeros(5,numreps);

% Add data to new matrix
QF(1,:) = qf0;

time = 2;
for i = 1:r-1

    curr_rep = data(i,1);

    if isnan(curr_rep)
        break
    end

    QF(time,curr_rep) = data(i,4);

    if time < 5 
        time = time+1;
    else
        time = 2;
    end

end

end
