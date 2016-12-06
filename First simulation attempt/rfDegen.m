function [symbols] = rfDegen(recwave,channel)

BW =  channel.bandwidth;
Time = channel.ts;
carrier = 2*pi*channel.carrier_freq;
points = channel.pointsPerWave;

number_of_symbols = length(recwave)/points;
x = linspace(0,Time,points);

for n = 1:number_of_symbols
    symbols(n) = 2*mean(recwave((n-1)*points+1:n*points).*cos(carrier*x))+...
                2*mean(recwave((n-1)*points+1:n*points).*sin(carrier*x))*1i;
    
    %    rfwave((n-1)*points+1:n*points) = real(symbols(n))*cos(carrier*x)+imag(symbols(n))*sin(carrier*x);
end
