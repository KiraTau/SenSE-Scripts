%create the folder
folderName = 'Data';

% Check if the folder exists. 'dir' specifies we are looking for a directory.
if ~exist(folderName, 'dir')
    mkdir(folderName);
    fprintf('Folder "%s" created.\n', folderName);
else
    fprintf('Folder "%s" already exists.\n', folderName);
end

%extract data from xlsx file and save it to mat
test_subjet = 8;
filename = strcat('Motion_NSF',num2str(test_subjet),'.xlsx');
sn = sheetnames(filename);%read all the sheets
MotionData = [];
start_day = [];
for j = 1:length(sn)
    currentName = sn(j);
    
    safeName = matlab.lang.makeValidName(currentName);
    
    % Load the data
    mtnData = readtable(filename, 'Sheet', currentName);
    
    %record the start day
    if isempty(start_day)
        start_day = datetime(mtnData{1,1});
    end

    %the columns are time in days, activity and light
    if size(mtnData,2) == 3%it does not have light recording, fill in 0
        MTN_Data = [table2array(mtnData(:,[2:3])) zeros(size(mtnData,1),1)];
    else%
        MTN_Data = table2array(mtnData(:,[2:4]));
    end

    for k = 1:size(mtnData,1)
        t = MTN_Data(k,1);
        MTN_Data(k,1) = t+days(datetime(mtnData{k,1}) - start_day);
    end
    %fill the time gap between with 0
    if ~isempty(MotionData)
        last_timestamp = MotionData(end, 1);
        current_first_timestamp = MTN_Data(1, 1);
        time_gap = current_first_timestamp - last_timestamp;
        minute_val = 1/(24*60*2);
        if time_gap > minute_val
            gap_times = (last_timestamp + minute_val : minute_val : current_first_timestamp - minute_val)';
            
            if ~isempty(gap_times)
                filler = [gap_times, zeros(length(gap_times), 2)];
                MotionData = [MotionData; filler];
            end
        end
    end
    MotionData = [MotionData;MTN_Data];
end

save(strcat('Data/',strcat('A',num2str(test_subjet)),'_MTN'),'MotionData','start_day')
