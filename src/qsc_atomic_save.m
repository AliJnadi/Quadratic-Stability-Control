function qsc_atomic_save(file, values)
% Local small checkpoints; MAT files contain only explicitly supplied values.
[folder,~,~] = fileparts(file);
if ~exist(folder,'dir'), mkdir(folder); end
temporary = [tempname(folder) '.mat'];
save(temporary,'-struct','values','-v7');
movefile(temporary,file,'f');
end
