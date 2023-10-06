% Sample data (replace this with your actual data)
% take only last noYears
function months_data_avg= reshapeCellToArrayAvg(data2,year)
noYears=1;
% data2=NPP;
data = data2((year-1)*365+1:year*365);%rand(10, 365);

% days_in_month = [31 29 31 30  31   30 31   31  30 31  30  31];
start_day =     [1  32 60 91  121 152 182 213 244 274 305 335];
end_day =       [31 59 90 120 151 181 212 243 273 304 334 365];

% Reshape the data into a 3D array
reshaped_data = reshape(data, noYears, 365, 1);

% Create a cell array to store each month's data
% months_data = cell(noYears, 12);

% Split the data into months
for j=1:noYears
    for month = 1:12
        months_data{j,month} = reshaped_data(:, start_day(month):end_day(month));
    end
end

% Find the mean value per month for each year -average over 2nd dimension
meanCell = cellfun(@(x) mean(x, 2), months_data, 'UniformOutput', false);

% convert cell to 9x12 array (1:3,:)=(3:6,:)=(6:9,:)
meanArray = cell2mat(meanCell);
% keep only top 3 row, as they are repeated below
months_data_avg=meanArray(1:noYears,:); 

%%
% figure(8)
% clf(8)
% plot(months_data_avg(1,:))
% 
% 
% axis tight
% xlabel('time (days)')
% ylabel('NPP')
 