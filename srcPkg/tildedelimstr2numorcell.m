function outarray = tildedelimstr2numorcell(instring)
% tildedelimstr2numorcell: Convert a tilde-delimited string to a double or cell 
% array.
% 
%   outarray = tildedelimstr2numorcell(instring)
%
%   instring Tilde-delimited string to be converted to an array.
% 
%   outarray Nx1 cell or double array containing the N tilde-delimited
%            elements of instring.  If all of the elements of instring can 
%            be converted to numbers then outarray is a double array; 
%            otherwise it is a cell array. 
% 
% -- Patrick J. Sutton
%    2005.07.10

%----- Make sure input is actually a string.
if (~ischar(instring))
    error('Input must be a string.')
end

%----- Initialize variables.
outarray = [];
cellout = 0;  %-- default output to double array.

%----- Parse input string.
j=0;  %-- number of elements parsed
while(~isempty(instring))
    %----- Get next element and convert to double.
    [out, instring] = strtok(instring,'~');
    j=j+1;
    if (cellout==0)
        %----- Test whether element can be converted to a double.
        if (~isempty(str2num(out)))
            outarray = [outarray, str2num(out)];
        else
            %----- Must switch output format to cell array.
            outarray = num2cell(outarray);
            cellout = 1;
        end
    end
    if (cellout==1)
        outarray{j} = out; 
    end
end

%----- Done.
return;
