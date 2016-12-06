function [recWave] = RFchannel(RFwave,channel)

switch(channel.type)
    case 'None'
        recWave = RFwave;
    case 'AWGN'
        recWave = RFwave+channel.Noise*randn(1,length(RFwave));
    otherwise
end