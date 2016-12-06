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
        Bits = zeros(length(symbols)*2,1);
        for n = 1:2:length(Bits)
            if (real(symbols(ceil(n/2)))>0)
                Bits(n) = 1;
            end
            if (imag(symbols(ceil(n/2)))>0)
                Bits(n+1) = 1;
            end
        end
    case '8PSK'
        Bits = zeros(length(symbols)*3,1);
        phase = angle(symbols);
        for n = 1:length(phase)
            if phase(n) < 0 
                phase(n)= phase(n)+2*pi; 
            end
            if (phase(n) > pi/8 && phase(n) < 3*pi/8) 
                Bits(3*(n-1)+1:3*(n-1)+3) = [1 1 0];
            elseif (phase(n) > 3*pi/8 && phase(n) < 5*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [0 1 0];
            elseif (phase(n) > 5*pi/8 && phase(n) < 7*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [0 1 1];
            elseif (phase(n) > 7*pi/8 && phase(n) < 9*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [0 0 1];
            elseif (phase(n) > 9*pi/8 && phase(n) < 11*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [0 0 0];
            elseif (phase(n) > 11*pi/8 && phase(n) < 13*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [1 0 0];
            elseif (phase(n) > 13*pi/8 && phase(n) < 15*pi/8)
                Bits(3*(n-1)+1:3*(n-1)+3) = [1 0 1];
            else
                Bits(3*(n-1)+1:3*(n-1)+3) = [1 1 1];
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

