function ctx = qsc_context(mask,grid)
% Cache symbolic derivations in memory; do not save symbolic models to disk.
if nargin < 1, mask = 15; end
if nargin < 2, grid = 2; end
persistent keys contexts
key = sprintf('%d_%d',grid,mask);
if isempty(keys), keys = {}; contexts = {}; end
index = find(strcmp(keys,key),1);
if ~isempty(index), ctx = contexts{index}; return; end
if grid == 2, fs = dec2bin(mask,4)-'0'; else, fs = ones(1,grid^2); end
ctx = make_context(fs,grid);
% Reference integration traces are unnecessary for compact provenance.
for field = {'z_traj','t'}
 if isfield(ctx.reference,field{1}), ctx.reference = rmfield(ctx.reference,field{1}); end
end
keys{end+1} = key; contexts{end+1} = ctx;
end
