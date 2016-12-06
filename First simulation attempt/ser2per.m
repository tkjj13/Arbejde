function [parallel_bits] = ser2per(serial_bits, number_of_lines)

if (mod(length(serial_bits),number_of_lines)==0)

    parallel_bits = reshape(serial_bits,number_of_lines,length(serial_bits)/number_of_lines)';

else
    error('Number of bits can not be parallelyse');
end