%% Input
% @debug: if debug is 1, do fprintf();
function DbgPrint(debug, varargin)
    if debug == 1
        fprintf(varargin{:});
        % fprintf(varargin{});          % This code is wrong, an error of invalid file identifier will be reported
    end
end
