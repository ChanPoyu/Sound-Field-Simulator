clear all 
close all

%% 8 x 20 soundfield reproduction 

sound_source = load('20ch_SquareSine_reproduction.mat');  %% after FIR conv
sound_source = struct2cell(sound_source);
sound_source = real(cell2mat(sound_source));
sound_source = sound_source .* 500;

sizeOfsoundSource = size(sound_source);

f_upper = 2000;
LPF = dsp.LowpassFilter;
LPF.SampleRate = 48000;
LPF.StopbandFrequency = f_upper * 1.25;
LPF.PassbandFrequency = f_upper;
[num,den] = tf(LPF);

LP_sound_source = zeros(20, length(num) + sizeOfsoundSource(2) - 1);

for n = 1 : 1 : 20
    LP_sound_source(n , :) = 1000 * conv(num, sound_source(n, :));
end

soundsource_x = 150;
soundsource_y = 200;


%% set speaker
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

%% control point
center_x = 300;
center_y = 200;

r = 20;

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
    round(center_x + r * cos(angle(1))), round(center_y + r * sin(angle(1)));
];


%% set color map
color_step = 50;
mymap = zeros(color_step * 2 + 1, 3);

for n = 1 : 1 : color_step
    mymap(n, :) = [n / color_step - 1 / color_step, n / color_step - 1 / color_step, 1];   
end
mymap(color_step + 1, :) = [1, 1, 1];

for n = 1 : 1 : color_step
    mymap(n + color_step  + 1, :) = [1, (color_step - n) / color_step, (color_step - n) / color_step];
end

%% set FDTD

% モデルの各定数
nx = 300;                   % X 方向の空間セル数 [pixels]
ny = 200;                   % Y 方向の空間セル数 [pixels]
dx = 10.0e-3;               % 空間刻み [m]
fs = 48000;
dt = 1/fs;                  % 時間刻み [s]
nmax = 15000;                % 計算ステップ数 [回]
savestep = 10;              % 保存間隔（ステップ数）


% 媒質の各定数
rho = 1.293;                % 媒質の密度ρ [kg/m^3]
kappa = 142.0e3;            % 媒質の体積弾性率κ [Pa]
c = 340;                    % speed of sound in air
CFL = c * dt / dx;

% PML
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

% initiate sound pressure
Vx = zeros(nx+2*W+1, ny+2*W);      % x方向粒子速度 [m/s]
Vy = zeros(nx+2*W,  ny+2*W+1);      % y方向粒子速度 [m/s]
P  = zeros(nx+2*W,  ny+2*W);      % 音圧 [Pa]
Px = zeros(nx+2*W,  ny+2*W); 
Py = zeros(nx+2*W,  ny+2*W); 

Precord = zeros(1 , nmax);


% animation frame
F(10000) = struct('cdata',[],'colormap',[]);
%% start simulation soundfield

for n = 1 : 1 : nmax
    
    
    Vx(2:end-1, :) = Vx(2:end-1, :) - (dt / rho) * Vx(2:end-1, :) .* Rx(2:end, :) - (dt / rho) * ((P(1: end-1, :) - P(2:end, :)) / dx);
    Vy(:, 2:end-1) = Vy(:, 2:end-1) - (dt / rho) * Vy(:, 2:end-1) .* Ry(:, 2:end) - (dt / rho) * ((P(:, 1: end-1) - P(:, 2:end)) / dx);
    
    
    Px(:, :) = Px(:, :) - dt * times(Px(:, :), Rx(:, :)) / rho - dt * ((kappa / dx) * (Vx(1:end-1, :) - Vx(2:end, :)));
    Py(:, :) = Py(:, :) - dt * times(Py(:, :), Ry(:, :)) / rho - dt * ((kappa / dx) * (Vy(:, 1:end-1) - Vy(:, 2:end)));
    
    
    P(:, :) = (Px(:, :) + Py(:, :));
    
%     % speaker rendering (Lowpassed)
    P(sp_position(1, 1), sp_position(1, 2)) = P(sp_position(1, 1), sp_position(1, 2)) + LP_sound_source(1, n);
    P(sp_position(2, 1), sp_position(2, 2)) = P(sp_position(2, 1), sp_position(2, 2)) + LP_sound_source(2, n);
    P(sp_position(3, 1), sp_position(3, 2)) = P(sp_position(3, 1), sp_position(3, 2)) + LP_sound_source(3, n);
    P(sp_position(4, 1), sp_position(4, 2)) = P(sp_position(4, 1), sp_position(4, 2)) + LP_sound_source(4, n);
    P(sp_position(5, 1), sp_position(5, 2)) = P(sp_position(5, 1), sp_position(5, 2)) + LP_sound_source(5, n);
    P(sp_position(6, 1), sp_position(6, 2)) = P(sp_position(6, 1), sp_position(6, 2)) + LP_sound_source(6, n);
    P(sp_position(7, 1), sp_position(7, 2)) = P(sp_position(7, 1), sp_position(7, 2)) + LP_sound_source(7, n);
    P(sp_position(8, 1), sp_position(8, 2)) = P(sp_position(8, 1), sp_position(8, 2)) + LP_sound_source(8, n);
    P(sp_position(9, 1), sp_position(9, 2)) = P(sp_position(9, 1), sp_position(9, 2)) + LP_sound_source(9, n);
    P(sp_position(10, 1), sp_position(10, 2)) = P(sp_position(10, 1), sp_position(10, 2)) + LP_sound_source(10, n);
    P(sp_position(11, 1), sp_position(11, 2)) = P(sp_position(11, 1), sp_position(11, 2)) + LP_sound_source(11, n);
    P(sp_position(12, 1), sp_position(12, 2)) = P(sp_position(12, 1), sp_position(12, 2)) + LP_sound_source(12, n);
    P(sp_position(13, 1), sp_position(13, 2)) = P(sp_position(13, 1), sp_position(13, 2)) + LP_sound_source(13, n);
    P(sp_position(14, 1), sp_position(14, 2)) = P(sp_position(14, 1), sp_position(14, 2)) + LP_sound_source(14, n);
    P(sp_position(15, 1), sp_position(15, 2)) = P(sp_position(15, 1), sp_position(15, 2)) + LP_sound_source(15, n);
    P(sp_position(16, 1), sp_position(16, 2)) = P(sp_position(16, 1), sp_position(16, 2)) + LP_sound_source(16, n);
    P(sp_position(17, 1), sp_position(17, 2)) = P(sp_position(17, 1), sp_position(17, 2)) + LP_sound_source(17, n);
    P(sp_position(18, 1), sp_position(18, 2)) = P(sp_position(18, 1), sp_position(18, 2)) + LP_sound_source(18, n);
    P(sp_position(19, 1), sp_position(19, 2)) = P(sp_position(19, 1), sp_position(19, 2)) + LP_sound_source(19, n);
    P(sp_position(20, 1), sp_position(20, 2)) = P(sp_position(20, 1), sp_position(20, 2)) + LP_sound_source(20, n);
    
    
    %% sound source original
%     P(sp_position(1, 1), sp_position(1, 2)) = P(sp_position(1, 1), sp_position(1, 2)) + sound_source(1, n);
%     P(sp_position(2, 1), sp_position(2, 2)) = P(sp_position(2, 1), sp_position(2, 2)) + sound_source(2, n);
%     P(sp_position(3, 1), sp_position(3, 2)) = P(sp_position(3, 1), sp_position(3, 2)) + sound_source(3, n);
%     P(sp_position(4, 1), sp_position(4, 2)) = P(sp_position(4, 1), sp_position(4, 2)) + sound_source(4, n);
%     P(sp_position(5, 1), sp_position(5, 2)) = P(sp_position(5, 1), sp_position(5, 2)) + sound_source(5, n);
%     P(sp_position(6, 1), sp_position(6, 2)) = P(sp_position(6, 1), sp_position(6, 2)) + sound_source(6, n);
%     P(sp_position(7, 1), sp_position(7, 2)) = P(sp_position(7, 1), sp_position(7, 2)) + sound_source(7, n);
%     P(sp_position(8, 1), sp_position(8, 2)) = P(sp_position(8, 1), sp_position(8, 2)) + sound_source(8, n);
%     P(sp_position(9, 1), sp_position(9, 2)) = P(sp_position(9, 1), sp_position(9, 2)) + sound_source(9, n);
%     P(sp_position(10, 1), sp_position(10, 2)) = P(sp_position(10, 1), sp_position(10, 2)) + sound_source(10, n);
%     P(sp_position(11, 1), sp_position(11, 2)) = P(sp_position(11, 1), sp_position(11, 2)) + sound_source(11, n);
%     P(sp_position(12, 1), sp_position(12, 2)) = P(sp_position(12, 1), sp_position(12, 2)) + sound_source(12, n);
%     P(sp_position(13, 1), sp_position(13, 2)) = P(sp_position(13, 1), sp_position(13, 2)) + sound_source(13, n);
%     P(sp_position(14, 1), sp_position(14, 2)) = P(sp_position(14, 1), sp_position(14, 2)) + sound_source(14, n);
%     P(sp_position(15, 1), sp_position(15, 2)) = P(sp_position(15, 1), sp_position(15, 2)) + sound_source(15, n);
%     P(sp_position(16, 1), sp_position(16, 2)) = P(sp_position(16, 1), sp_position(16, 2)) + sound_source(16, n);
%     P(sp_position(17, 1), sp_position(17, 2)) = P(sp_position(17, 1), sp_position(17, 2)) + sound_source(17, n);
%     P(sp_position(18, 1), sp_position(18, 2)) = P(sp_position(18, 1), sp_position(18, 2)) + sound_source(18, n);
%     P(sp_position(19, 1), sp_position(19, 2)) = P(sp_position(19, 1), sp_position(19, 2)) + sound_source(19, n);
%     P(sp_position(20, 1), sp_position(20, 2)) = P(sp_position(20, 1), sp_position(20, 2)) + sound_source(20, n);
 
    Precord(n) = P(mic_position(7, 1), mic_position(7, 2));
   
    
    if(n > 3500 && n <= 3500 + 10000)
        % visualize
        clims = [-5 5];
        imagesc((P(101:400, 101: 300))', clims);     % 行列の画像表示（正負の値を表示する場合）：↑いずれかのみ
        
        axis equal; axis tight;             % 軸の調整
        title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
        xlabel('x'); ylabel('y');
        colormap(mymap);
        colorbar;
        hold on;
        scatter(sp_position(:, 1) - 100, sp_position(:, 2) - 100, 'filled', 'k');
        plotmic = plot(mic_position(:, 1) - 100, mic_position(:, 2) - 100, 'k');
        plotmic.Marker = '*';
        plotmic.LineWidth = 2;
%         source_pos = plot(soundsource_x - 100, soundsource_y - 100, 'k');
%         source_pos.Marker = 'diamond';

%         drawnow limitrate;
        F(n - 3500) = getframe(gcf);
        if(rem(n, 80)  == 0)
            close all
        end
        
    end
    
    
   
end

movie(F,2);


plot(Precord);
