function [bin] = myDec2bin(dec)

s = length(dec);
rows = 1;
dec_max = max(dec);
while (dec_max > 1)
    rows = rows+1;
    dec_max = (dec_max-mod(dec_max,2))/2;
end

    bin = zeros(s(1),rows);
for n = 1:s(1)
    for k = rows:-1:1
        bin(n,k) = mod(dec(n),2);
        dec(n) = (dec(n)-bin(n,k))/2;
    end
end