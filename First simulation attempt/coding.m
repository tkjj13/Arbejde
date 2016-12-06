function [code_bits] = coding(Bits, crypto_scheme, soruce_scheme, channel_scheme)


switch(crypto_scheme)
    case 'None'
        crypt_bits = Bits;
    case 'User_defined'
        crypt_bits = Bits;
    otherwise     
        error('No crypto coding scheme of that name exits')
end

switch(soruce_scheme)
    case 'None'
        source_bits = crypt_bits;
    case 'User_defined'
        source_bits = crypt_bits;
    otherwise     
        error('No source coding scheme of that name exits')
end


switch(channel_scheme)
    case 'None'
        code_bits = source_bits;
    case 'User_defined'
        code_bits = source_bits;
    otherwise     
        error('No channel coding scheme of that name exits')
end


