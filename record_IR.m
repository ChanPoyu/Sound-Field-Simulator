%% FDTD PML

clear all;
close all;

% モデルの各定数
nx = 300;                   % X 方向の空間セル数 [pixels]
ny = 200;                   % Y 方向の空間セル数 [pixels]
dx = 10.0e-3;               % 空間刻み [m]
fs = 48000;
dt = 1/fs;                  % 時間刻み [s]
nmax = 10000;                % 計算ステップ数 [回]
savestep = 10;              % 保存間隔（ステップ数）
nIinterval = 100;

% 媒質の各定数
rho = 1.293;                % 媒質の密度ρ [kg/m^3]
kappa = 142.0e3;            % 媒質の体積弾性率κ [Pa]
c = 340;                    % speed of sound in air
CFL = c * dt / dx;

% PML Layers
W = 100; 

Rx = zeros(nx+2*W,  ny+2*W); 
Ry = zeros(nx+2*W,  ny+2*W);
Rmax = 10 * rho * c;
m = 6;

for n = 1 : 1 : W
   Rx(n, :) = Rmax * ((W+1 - n) / W)^m;
   Rx(nx+2*W+1-n, :) = Rx(n, :);
end

for n = 1 : 1 : W
   Ry(:, n) = Rmax * ((W+1 - n) / W)^m;
   Ry(:, ny+2*W + 1 - n) = Ry(:, n);
end


% 点音源の各定数
sig_freq = 100;           % 印可波形の周波数 [Hz]
sig_amp = 3000.0;           % 印可波形の振幅
sig_x = (150);                % 点音源の x 座標 [pixels] ( 0 < sig_x <= nx)
sig_y = (ny + 2 * W) / 2;                % 点音源の y 座標 [pixels] ( 0 < sig_y <= ny)


% initiate sound pressure
Vx = zeros(nx+2*W+1, ny+2*W);      % x方向粒子速度 [m/s]
Vy = zeros(nx+2*W,  ny+2*W+1);      % y方向粒子速度 [m/s]
P  = zeros(nx+2*W,  ny+2*W);      % 音圧 [Pa]
Px = zeros(nx+2*W,  ny+2*W); 
Py = zeros(nx+2*W,  ny+2*W); 



%% load TSP

y = load("arbitrary_TSP.mat");
y = struct2cell(y);
y = 300 * real(cell2mat(y));


%% Gaussian pulse

% a = 3;
% for xx = 1 : 1 : nx+2*W
%     for yy = 1 : 1 : ny+2*W
%         r_sqr = (xx - sig_x)^2 + (yy - sig_y)^2;
%         if r_sqr < 10000
%             P(xx, yy) = sig_amp * exp(- r_sqr / a^2);
%         end
%     end
% end


%% position of speakers

x_lin = linspace(100, 400, 7);
y_lin = linspace(100, 300, 5);

sp_position = [
    x_lin(1), y_lin(1);
    x_lin(2), y_lin(1);
    x_lin(3), y_lin(1);
    x_lin(4), y_lin(1);
    x_lin(5), y_lin(1);
    x_lin(6), y_lin(1);
    x_lin(7), y_lin(1);
    x_lin(7), y_lin(2);
    x_lin(7), y_lin(3);
    x_lin(7), y_lin(4);
    x_lin(7), y_lin(5);
    x_lin(6), y_lin(5);
    x_lin(5), y_lin(5);
    x_lin(4), y_lin(5);
    x_lin(3), y_lin(5);
    x_lin(2), y_lin(5);
    x_lin(1), y_lin(5);
    x_lin(1), y_lin(4);
    x_lin(1), y_lin(3);
    x_lin(1), y_lin(2);
];

%% position of mic

center_x = 300;
center_y = 200;

r = 20;

dd = 2 * pi * r * dx / 8;

f_limit = c / (2 * dd);

angle = [pi / 2, pi / 4, 0, -pi / 4, -pi / 2, -pi * 3 / 4, -pi, pi * 3 / 4];

mic_position = [
    round(center_x + r * cos(angle(1))), round(center_y + r * sin(angle(1)));
    round(center_x + r * cos(angle(2))), round(center_y + r * sin(angle(2)));
    round(center_x + r * cos(angle(3))), round(center_y + r * sin(angle(3)));
    round(center_x + r * cos(angle(4))), round(center_y + r * sin(angle(4)));
    round(center_x + r * cos(angle(5))), round(center_y + r * sin(angle(5)));
    round(center_x + r * cos(angle(6))), round(center_y + r * sin(angle(6)));
    round(center_x + r * cos(angle(7))), round(center_y + r * sin(angle(7)));
    round(center_x + r * cos(angle(8))), round(center_y + r * sin(angle(8)));
];



Precord = zeros(8,20,nmax);

mic_record = zeros(8, nmax);


%% color map configure
color_step = 50;
mymap = zeros(color_step * 2 + 1, 3);

for n = 1 : 1 : color_step
    mymap(n, :) = [n / color_step - 1 / color_step, n / color_step - 1 / color_step, 1];   
end
mymap(color_step + 1, :) = [1, 1, 1];

for n = 1 : 1 : color_step
    mymap(n + color_step  + 1, :) = [1, (color_step - n) / color_step, (color_step - n) / color_step];
end

start = cputime;

%% start calculate 
for speaker_index = 1 : 1 : 20



    for n = 1 : 1 : nmax
        
        
        Vx(2:end-1, :) = Vx(2:end-1, :) - (dt / rho) * Vx(2:end-1, :) .* Rx(2:end, :) - (dt / rho) * ((P(1: end-1, :) - P(2:end, :)) / dx);
        Vy(:, 2:end-1) = Vy(:, 2:end-1) - (dt / rho) * Vy(:, 2:end-1) .* Ry(:, 2:end) - (dt / rho) * ((P(:, 1: end-1) - P(:, 2:end)) / dx);
        
        
        Px(:, :) = Px(:, :) - dt * times(Px(:, :), Rx(:, :)) / rho - dt * ((kappa / dx) * (Vx(1:end-1, :) - Vx(2:end, :)));
        Py(:, :) = Py(:, :) - dt * times(Py(:, :), Ry(:, :)) / rho - dt * ((kappa / dx) * (Vy(:, 1:end-1) - Vy(:, 2:end)));
        
        
        P(:, :) = (Px(:, :) + Py(:, :));
        
        
        if(n <= length(y))
            P(sp_position(speaker_index, 1), sp_position(speaker_index, 2)) = P(sp_position(speaker_index, 1), sp_position(speaker_index, 2)) + y(n);
        end
        
        Precord(speaker_index, 1, n) = P(mic_position(1, 1), mic_position(1, 2));
        Precord(speaker_index, 2, n) = P(mic_position(2, 1), mic_position(2, 2));
        Precord(speaker_index, 3, n) = P(mic_position(3, 1), mic_position(3, 2));
        Precord(speaker_index, 4, n) = P(mic_position(4, 1), mic_position(4, 2));
        Precord(speaker_index, 5, n) = P(mic_position(5, 1), mic_position(5, 2));
        Precord(speaker_index, 6, n) = P(mic_position(6, 1), mic_position(6, 2));
        Precord(speaker_index, 7, n) = P(mic_position(7, 1), mic_position(7, 2));
        Precord(speaker_index, 8, n) = P(mic_position(8, 1), mic_position(8, 2));
        
        if rem(n, nIinterval) == 0
            display(speaker_index);
            display(n);
        end
        
%         if n < (1/sig_freq)/dt
%             sig = sig_amp * sin(2.0*pi*sig_freq*n*dt) ^ 5;
% 
%             P(sig_x,sig_y) = P(sig_x,sig_y) + sig;
%         end

        % record sound source
%         mic_record(1, n) = P(mic_position(1, 1), mic_position(1, 2));
%         mic_record(2, n) = P(mic_position(2, 1), mic_position(2, 2));
%         mic_record(3, n) = P(mic_position(3, 1), mic_position(3, 2));
%         mic_record(4, n) = P(mic_position(4, 1), mic_position(4, 2));
%         mic_record(5, n) = P(mic_position(5, 1), mic_position(5, 2));
%         mic_record(6, n) = P(mic_position(6, 1), mic_position(6, 2));
%         mic_record(7, n) = P(mic_position(7, 1), mic_position(7, 2));
%         mic_record(8, n) = P(mic_position(8, 1), mic_position(8, 2));
%         
        
% if(n > 0)
%     clims = [-5 5];
%     imagesc((P)', clims);     % 行列の画像表示（正負の値を表示する場合）：↑いずれかのみ
%     
%     axis equal; axis tight;             % 軸の調整
%     title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
%     xlabel('x'); ylabel('y');
%     colormap(mymap);
%     colorbar;
%     hold on;
%     scatter(sp_position(:, 1), sp_position(:, 2), 'filled', 'k');
%     plotmic = plot(mic_position(:, 1), mic_position(:, 2), 'k');
%     plotmic.Marker = '*';
%     plotmic.LineWidth = 2;
%     drawnow;
% end
%         
        
        
        
    end



end

endTime = cputime - start;
display(endTime);

iTSP = load("arbitrary_TSP_inverse.mat");
iTSP = struct2cell(iTSP);
iTSP = real(cell2mat(iTSP));

IR = zeros(20, 8, nmax + length(iTSP) - 1);


for mic_index = 1 : 1 : 8
    for sp_index = 1 : 1 : 20
        
        TSP_record = reshape(Precord(sp_index, mic_index, :), [1, nmax]);
        
        IR(sp_index, mic_index, :) = conv(TSP_record, iTSP);

    end
end


figure(1);
subplot(1,2,1);
plot(reshape(Precord(1,1,:), [1, nmax]));
subplot(1,2,2);
plot(reshape(IR(1,1,:), [1, nmax + length(iTSP) - 1]));


save('IR_freeSoundField.mat', 'IR');

% plot(mic_record(7, :));
% save('8ch_SineSquare_record.mat', 'mic_record');




