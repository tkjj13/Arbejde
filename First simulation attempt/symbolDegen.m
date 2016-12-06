function [Bits] = symbolDegen(symbols, scheme)

switch(scheme)
    case 'BPSK'
        Bits = zeros(length(symbols),1);
        for n = 1:length(symbols)
            if (symbols(n)>0)
                Bits(n) = 1;
            end
        end
    case 'QPSK'
        Bits = zeros(length(symbols)/2,2);
        for n = 1:length(Bits)
            if (real(symbols(n))>0)
                Bits(n,1) = 1;
            end
            if (imag(symbols(n))>0)
                Bits(n,2) = 1;
            end
        end
    case '8PSK'
        parallel_bits = ser2per(Bits,3);
        [number_of_symbols, ~] = size(parallel_bits);
        symbols = zeros(number_of_symbols,1);
        for n = 1:number_of_symbols
            switch(bin2dec(sprintf('%i%i%i',parallel_bits(n,1),parallel_bits(n,2),parallel_bits(n,3))))
                case 7
                    symbols(n) = exp(1i*0/8*2*pi);
                case 6
                    symbols(n) = exp(1i*1/8*2*pi);
                case 2
                    symbols(n) = exp(1i*2/8*2*pi);
                case 3
                    symbols(n) = exp(1i*3/8*2*pi);
                case 1
                    symbols(n) = exp(1i*4/8*2*pi);
                case 0
                    symbols(n) = exp(1i*5/8*2*pi);
                case 4
                    symbols(n) = exp(1i*6/8*2*pi);
                case 5
                    symbols(n) = exp(1i*7/8*2*pi);
                otherwise
            end
        end
    case '64QAM'
        parallel_bits = ser2per(Bits,6);
        [number_of_symbols, ~] = size(parallel_bits);
        symbols = zeros(number_of_symbols,1);
        for n = 1:number_of_symbols
            switch(bin2dec(sprintf('%i%i%i',parallel_bits(n,4),parallel_bits(n,5),parallel_bits(n,6))))
                case 7
                    symbols(n) = 3;
                case 6
                    symbols(n) = 1;
                case 2
                    symbols(n) = -1;
                case 3
                    symbols(n) = -3;
                case 1
                    symbols(n) = -5;
                case 0
                    symbols(n) = -7;
                case 4
                    symbols(n) = 7;
                case 5
                    symbols(n) = 5;
                otherwise
            end
            switch(bin2dec(sprintf('%i%i%i',parallel_bits(n,1),parallel_bits(n,2),parallel_bits(n,3))))
                case 7
                    symbols(n) = symbols(n)-3*1i;
                case 6
                    symbols(n) = symbols(n)-1*1i;
                case 2
                    symbols(n) = symbols(n)+1*1i;
                case 3
                    symbols(n) = symbols(n)+3*1i;
                case 1
                    symbols(n) = symbols(n)+5*1i;
                case 0
                    symbols(n) = symbols(n)+7*1i;
                case 4
                    symbols(n) = symbols(n)-7*1i;
                case 5
                    symbols(n) = symbols(n)-5*1i;
                otherwise
            end
        end
    otherwise
        error('Not a valid modulation scheme');
end

