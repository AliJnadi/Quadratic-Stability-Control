function manifest = qsc_manifest(root,cfg)
files = [dir(fullfile(root,'*.m')); dir(fullfile(root,'src','*.m')); ...
 dir(fullfile(root,'experiments','*.m')); dir(fullfile(root,'analysis','*.m'))];
records = struct('path',{},'sha256',{});
for j = 1:numel(files)
 file = fullfile(files(j).folder,files(j).name);
 fid = fopen(file,'rb'); bytes = fread(fid,Inf,'*uint8'); fclose(fid);
 records(j).path = file(numel(root)+2:end);
 records(j).sha256 = hash_bytes(bytes);
end
[~,order] = sort({records.path}); records = records(order);
% cvx_version reloads preferences and can discard initialized solver callbacks.
% Preserve the live CVX state while collecting version metadata.
global cvx___
cvx_state = cvx___;
try
 cvx_text = evalc('cvx_version');
catch exception
 cvx___ = cvx_state;
 rethrow(exception);
end
cvx___ = cvx_state;
cvx_version_line = regexp(cvx_text,'Version [^\r\n]+','match','once');
environment = struct('matlab',version,'computer',computer, ...
 'cvx',cvx_version_line,'toolboxes',ver);
manifest = struct('schema',1,'created',datestr(now,30),'configuration',cfg, ...
 'environment',environment,'sources',records);
manifest.signature = hash_bytes(unicode2native(jsonencode( ...
 struct('configuration',cfg,'environment',environment,'sources',records)),'UTF-8'));
end
function h = hash_bytes(bytes)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(uint8(bytes(:)),'int8'));
h = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
