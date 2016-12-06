function [dec] = myBin2dec(bin)

s = size(bin);
dec = zeros(s(1),1);
for n = 1:s(1)
    for k = 1:s(2)
        dec(n) = dec(n)+2^(s(2)-k)*bin(n,k);
    end
end