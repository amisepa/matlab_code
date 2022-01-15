%%This function finds the egg file's match from IDL database
function match = findMatchIDL(id_eeg, list_ids_idl, date_eeg, idl_startDate, idl_endDate)

% get all IDs to same format
for a = 1:size(list_ids_idl,1)
    id = char(list_ids_idl(a));
    if length(id) < 8
        id(:,end+1:10) = 'x';
    elseif length(id) > 10
        id = id(:,1:10);
    end
    list_ids(a,:) = {id};
end

matchID = contains(list_ids,id_eeg);

if sum(matchID) > 1
    count = 1;
    for candidates = 1:size(list_ids,1)
        if true(matchID(candidates))
            match(count) = candidates;
            if ~isnat(idl_endDate(candidates))
                date_idl(count,:) = datetime(idl_endDate(candidates),'Format','y-MMM-d HH:mm:ss');
                count = count+1;
            else
                date_idl(count,:) = datetime(idl_startDate(candidates),'Format','y-MMM-d HH:mm:ss');
                count = count+1;
            end
        end
    end
    
    %Find which date is closest to find the right match
    [~,ind] = min(abs(datenum(date_idl)-datenum(date_eeg)));
    match = match(:,ind);
    
elseif sum(matchID) == 1
    match = find(matchID == 1); 
    
else
    match = NaN;
    
end

