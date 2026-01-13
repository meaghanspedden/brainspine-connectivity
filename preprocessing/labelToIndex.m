function idxCell = labelToIndex(label)
    % label: string like 'A2' or a cell array {'C1','D3',...}

    % Ensure we can index into it uniformly
    if ischar(label) || isstring(label)
        label = {char(label)};
    end

    n = numel(label);
    idxStr = strings(1,n);

    for k = 1:n
        this = label{k};
        colChar = upper(this(1));          % A–H
        row     = str2double(this(2:end)); % 1–8 (robust if "A10" ever appears)

        col = colChar - 'A' + 1;           % convert A–H → 1–8

        % Compute linear index (column-major)
        idx = (col - 1) * 8 + row;

        % Store cell array
         idxCell{k} = num2str(idx);
    end
end