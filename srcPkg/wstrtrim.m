function s = wstrtrim(s)
% remove leading and trailing whitespace from a string

if (~isempty(s))
    while (isspace(s(1)))
        s = s(2:end);
    end

    while (isspace(s(end)))
        s = s(1:(end - 1));
    end
end