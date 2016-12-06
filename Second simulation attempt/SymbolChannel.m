function [rec_symbols] = SymbolChannel(symbols,channel)

switch(channel.type)
    case 'Reyleigh'
        rec_symbols = symbols;
    case 'AWGN'
        rec_symbols = symbols+channel.Noise*randn(length(symbols),1)+channel.Noise*randn(length(symbols),1)*1i;
    otherwise
end