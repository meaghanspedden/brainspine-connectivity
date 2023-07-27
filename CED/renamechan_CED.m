function chantitle  = renamechan_CED(chantitle,channr,otherchanstitles,printinfos)

% renamechan_CED.m est utilisée par la fonction readCEDevents.m :
% - elle remplace les caractères spéciaux, interdits dans les noms de variables,
% par des underscores '_'.
% - si le canal n'a pas de nom, elle le remplace par le numéro du canal
% Spike2
% - si plusieurs canaux ont le même nom, ils sont différenciés en ajoutant le
% numéro de canal Spike2
% 

% cant make fieldnames with minusses
if any(ismember('-',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with minusses)\n', chantitle, strrep(chantitle,'-','_'));
    end
    chantitle = strrep(chantitle,'-','_');
end

% cant make fieldnames with points
if any(ismember('.',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with points)\n', chantitle, strrep(chantitle,'.','_'));
    end
    chantitle = strrep(chantitle,'.','_');
end

% cant make fieldnames with parentheses
if any(ismember('(',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with parentheses)\n', chantitle, strrep(chantitle,'(','_'));
    end
    chantitle = strrep(chantitle,'(','_');
end

% cant make fieldnames with parentheses
if any(ismember(')',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with parentheses)\n', chantitle, strrep(chantitle,')','_'));
    end
    chantitle = strrep(chantitle,')','_');
end

% cant make fieldnames with white spaces
if any(ismember(' ',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with white spaces)\n', chantitle, strrep(chantitle,' ','_'));
    end
    chantitle = strrep(chantitle,' ','_');
end

% cant make fieldnames with white /
if any(ismember('/',chantitle))
    if printinfos
        fprintf('Channel %s is renamed %s (cant make fieldnames with '/')\n', chantitle, strrep(chantitle,'/','_'));
    end
    chantitle = strrep(chantitle,'/','_');
end

%create name if is empty
if isempty(chantitle)
    chantitle = sprintf('chan%d', channr);
    if printinfos
        fprintf('Channel %d is renamed %s (has no name)\n', channr, chantitle);
    end
end

%add 'x' to the field name if it begins with a number or a special
%character (replaced by '_' above)
if ismember(chantitle(1), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '_'])
    chantitle = insertBefore(chantitle,chantitle,'X');
    if printinfos
        fprintf('Channel %s is renamed %s (begins with a number)\n', chantitle(2:end), chantitle);
    end
end

%rename channel if it has the same name that a previous one
if any(strcmp(chantitle, otherchanstitles))
    oldname = chantitle;
    chantitle = sprintf('%s_chan%d', chantitle, channr);
    if printinfos
        fprintf('Channel %s (%d) as the same name than a previous channel : renamed %s\n', oldname, channr, chantitle);
    end
end

end