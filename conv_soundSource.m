clear all
close all

%% generate sound source for sound field reproduction;


% [sound_source, fs] = audioread();

% FIR Filter M x N channels
Inverse_Filter = load('INVERSE_FILTER.mat');
Inverse_Filter = struct2cell(Inverse_Filter);
Inverse_Filter = real(cell2mat(Inverse_Filter));
sizeOfInversefilter = size(Inverse_Filter);


% load sound source 
Sound_source = load('8ch_SineSquare_record.mat');  % mic signal M channels
Sound_source = struct2cell(Sound_source);
Sound_source = real(cell2mat(Sound_source));
sizeOfSoundsource = size(Sound_source);
% plot(Sound_source(7, :));
% hold on;
% plot(Sound_source(1, :));
% plot(Sound_source(2, :));
% plot(Sound_source(3, :));

% output signal
Output_signal = zeros(20, sizeOfSoundsource(2) + sizeOfInversefilter(3) - 1);
sizeofOutputSignal = size(Output_signal);

for ii = 1 : 1 : sizeofOutputSignal(1)
    for jj = 1 : 1 : sizeOfSoundsource(1)
        FIR = reshape(Inverse_Filter(jj, ii, :), [sizeOfInversefilter(3), 1]);
        Output_signal(ii, :) = Output_signal(ii, :) + conv(FIR, Sound_source(jj, :));
    end
end

plot(Output_signal(1 , :));
save('20ch_SquareSine_reproduction.mat', 'Output_signal');