function qsc_write_json(file, value)
fid = fopen(file,'w'); assert(fid >= 0, 'Cannot write %s',file);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end
