function s = basename(string)
% BASENAME return base name of path
%    S = BASENAME(STRING) returns everything after the last '/' in STRING

s = regexprep(string, '.*/', '');
