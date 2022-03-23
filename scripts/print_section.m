function [] = print_section(title, n_chars)
%FUNCTION [] = print_section(title)
% Prints a framed section header into console

if nargin < 1
    title = 'New Section';
end
if nargin < 2
    n_chars = 75;
end

% add datetime
title = char([datestr(now) ' ' title]);

% adjust length
if length(title) > n_chars-3
    n_chars = length(title)+3;
end

% get number of spaces
lineend = repmat(' ', 1, n_chars-length(title)-3);

% print
fprintf('\n┌%s┐\n┝ %s%s┤\n└%s┘\n\n', ...
        repmat('-', 1, n_chars-2), title, lineend, repmat('-', 1, n_chars-2))

end
