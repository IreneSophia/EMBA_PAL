clearvars
data = struct();

% get a list of preprocessed csv data files
files = dir("PAL*data.csv");

% loop through them
for i = 1:length(files)
    % read in the data
    prepro = readtable(files(i).name);
    % recode difficulty to noise
    prepro.noise(strcmp(prepro.difficulty, "easy"))      = 0.1;
    prepro.noise(strcmp(prepro.difficulty, "medium"))    = 0.3;
    prepro.noise(strcmp(prepro.difficulty, "difficult")) = 0.9;
    % get all the subIDs
    subIDs = unique(prepro.subID);
    % loop through subjects
    for s = 1:length(subIDs)
        % get this subs data
        sub = prepro(prepro.subID == subIDs(s),:);
        if height(sub) ~= 336
            warning('Subject %d has wrong trial number: %d', subIDs(s), height(sub))
        end
        % extract ut
        data(s).u = [sub.ut sub.noise];
        % exctract correct, log transformed rts
        sub.rt(strcmp(sub.acc, 'FALSE')) = nan(sum(strcmp(sub.acc, 'FALSE')), 1);
        data(s).y = log(sub.rt);
        % extract group
        data(s).group = sub.diagnosis{1};
        % add subID
        data(s).subID = subIDs(s);
    end
    
    save(strrep(files(i).name, 'csv', 'mat'), 'data')

end