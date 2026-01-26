function vprint(verbose, level, varargin)
    if verbose >= level
        fprintf(varargin{:});
    end
end