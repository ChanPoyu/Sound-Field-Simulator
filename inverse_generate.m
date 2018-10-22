clear all
close all

%% Read IR and FFT


ir = load("IR_freeSoundField.mat");
ir = struct2cell(ir);
ir = real(cell2mat(ir));





%% cut IR 5000 ~ 13192

ir_cut = ir(:, : , 5000 : 13192);

sizeOfIr = size(ir_cut);

IR = zeros(sizeOfIr);

for micIndex = 1 : 1 : sizeOfIr(1)
    for speakerIndex = 1 : 1 : sizeOfIr(2)
        signal = reshape( ir_cut(micIndex, speakerIndex, :), [1, sizeOfIr(3)]);
        SIGNAL = fft(signal);
        IR(micIndex, speakerIndex, :) = SIGNAL;
    end
end

% semilogx(abs(reshape( IR(1, 1, :), [1, sizeOfIr(3)])));

%% Inverse

H = zeros(sizeOfIr(2), sizeOfIr(1), sizeOfIr(3));
lambda = 1;
I = eye(sizeOfIr(2))';

for FrequencyBin_Index = 1 : 1 : sizeOfIr(3)
    G = IR(:, :, FrequencyBin_Index);
%     H(:, :, FrequencyBin_Index) = pinv(G);

    H(:, :, FrequencyBin_Index) = (G' * G + lambda * I) \ I * G';


end

%% make FIR

bin = 257;

FIR_filter = zeros(sizeOfIr(2), sizeOfIr(1), sizeOfIr(3));

for  speakerIndex= 1 : 1 : sizeOfIr(2)
    for micIndex = 1 : 1 : sizeOfIr(1)
        freq_bin = reshape(H(speakerIndex, micIndex, :), [1, sizeOfIr(3)]);
        FIR_filter(speakerIndex, micIndex, :) = circshift(real(ifft(freq_bin)), round(sizeOfIr(3) / 2) );
    end
end


% centerOfFIR = round(sizeOfIr(3) / 2) ;
% FIR_filter = FIR_filter(: , : , centerOfFIR - 2000 : centerOfFIR + 2000);


% hamming_window = hamming(4001);
% 
% for  speakerIndex= 1 : 1 : sizeOfIr(2)
%     for micIndex = 1 : 1 : sizeOfIr(1)
%         FIR_temp = reshape(FIR_filter(speakerIndex, micIndex, :), [4001 1]);
%         FIR_temp = FIR_temp .* hamming_window;
%         FIR_filter(speakerIndex, micIndex, :) = FIR_temp;
%     end
% end



plot(reshape(FIR_filter(1, 4, :), [1, 8193]));

save('INVERSE_FILTER.mat', 'FIR_filter');
    