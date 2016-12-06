clc;
close all;
clear all;


schemes = struct;

schemes.crypto = 'None';
schemes.source = 'None';
schemes.channel = 'None';
schemes.mod = '64QAM';
channel = struct;    
channel.type = 'AWGN';
channel.Noise = 0.3;
channel.bandwidth = 200E3;
channel.stopBandAtt = 30;
channel.ts = 1E-4;
channel.carrier_freq = 2E6;
channel.pointsPerWave = 100;

users = 2;     % number of user in the system
room_size = 1;  % length of square room
n = 9996;
Bits = round(rand(n,1));
[coded_bits] = coding(Bits, schemes.crypto, schemes.source, schemes.channel);
[symbols] = symbolGen(coded_bits, schemes.mod);
[rfwave] = rfgen(symbols,channel);
[recwave] = RFchannel(rfwave, channel);
%[rec_symbols] = SymbolChannel(symbols, channel);
[rec_symbols] = rfDegen(recwave,channel);
figure;
plot(symbols,'*')
hold on
plot(rec_symbols,'*')
grid
figure;
plot(rfwave)
return
% 
BS_placement = [room_size/2 room_size/2];
user_placement = unifrnd(0,room_size,users,2);


for n = 1:users
    direct_segments(n,:) = [BS_placement(1) user_placement(n,1) BS_placement(2) user_placement(n,2)];
end

% Reflections

for n = 1:users
    ref_point_left = [0 (BS_placement(2)-user_placement(n,2))/((BS_placement(1)/user_placement(n,1))+1)+user_placement(n,2)];
    ref_point_right = [room_size (BS_placement(2)-user_placement(n,2))/(((room_size-BS_placement(1))/(room_size-user_placement(n,1)))+1)+user_placement(n,2)];
    ref_point_top = [(BS_placement(1)-user_placement(n,1))/(((room_size-BS_placement(2))/(room_size-user_placement(n,2)))+1)+user_placement(n,1) room_size];
    ref_point_bottom = [(BS_placement(1)-user_placement(n,1))/((BS_placement(2)/user_placement(n,2))+1)+user_placement(n,1) 0];
    
    left_segments(n,:) = [BS_placement(1) ref_point_left(1) user_placement(n,1) BS_placement(2) ref_point_left(2) user_placement(n,2)];
    right_segments(n,:) = [BS_placement(1) ref_point_right(1) user_placement(n,1) BS_placement(2) ref_point_right(2) user_placement(n,2)];
    top_segments(n,:) = [BS_placement(1) ref_point_top(1) user_placement(n,1) BS_placement(2) ref_point_top(2) user_placement(n,2)];
    bottom_segments(n,:) = [BS_placement(1) ref_point_bottom(1) user_placement(n,1) BS_placement(2) ref_point_bottom(2) user_placement(n,2)];
end


scatter(user_placement(:,1),user_placement(:,2));
hold on;
scatter(BS_placement(1),BS_placement(2),'r')
for n = 1:users
    plot(direct_segments(n,1:2),direct_segments(n,3:4),'r');
    plot(left_segments(n,1:3),left_segments(n,4:6),'g');
    plot(right_segments(n,1:3),right_segments(n,4:6),'g');
    plot(top_segments(n,1:3),top_segments(n,4:6),'g');
    plot(bottom_segments(n,1:3),bottom_segments(n,4:6),'g');
end
    xlim([0 room_size]);
ylim([0 room_size]);
grid;