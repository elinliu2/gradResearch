f = @() timing(); % handle to function
timeit(f)

function [] = timing()
    for i = 1:262143
        A = rand(41);
        decimalNumber = i;
        binaryStr = dec2bin(decimalNumber, 41);          % returns '1101' as a string
        b = binaryStr - '0';  
        x = A\b';
    end
end