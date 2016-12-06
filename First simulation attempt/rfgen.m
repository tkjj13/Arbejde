function [rfwave] = rfgen(symbols,channel)

number_of_symbols = length(symbols);
f0 =  channel.bandwidth/2;
Att = channel.stopBandAtt;
Time = channel.ts;
carrier = 2*pi*channel.carrier_freq;
points = channel.pointsPerWave;
%order = ceil(30/(20*log10(((BW/2)*sqrt(2))/((BW/2)/sqrt(2)))));
%order = ceil(Att/(20*log10(2))); % det samme som ovenover

% filter symbols 
   d = fdesign.lowpass('Fp,Fst,Ap,Ast',f0/sqrt(2),f0*sqrt(2),0.5,Att,points/Time);
   Hd = design(d,'equiripple');
   
% interpolate
   interSymbols = reshape(repmat(symbols',points,1),symbols*points,1);
   output = filter(Hd,real(symbols))+filter(Hd,imag(symbols))*1i;




x = linspace(0,Time,points);
for n = 1:number_of_symbols
    rfwave((n-1)*points+1:n*points) = ...
    real(symbols(n))*cos(carrier*x)+imag(symbols(n))*sin(carrier*x);
end
