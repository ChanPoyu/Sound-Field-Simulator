clear all;
close all;

%% FDTD PML

% モデルの各定数
nx = 300;                   % X 方向の空間セル数 [pixels]
ny = 200;                   % Y 方向の空間セル数 [pixels]
dx = 10.0e-3;               % 空間刻み [m]
fs = 48000;
dt = 1/fs;                  % 時間刻み [s]
nmax = 10000;                % 計算ステップ数 [回]
savestep = 10;              % 保存間隔（ステップ数）


% 媒質の各定数
rho = 1.293;                % 媒質の密度ρ [kg/m^3]
kappa = 142.0e3;            % 媒質の体積弾性率κ [Pa]
c = 340;                    % speed of sound in air
CFL = c * dt / dx;

R_wall = 0.8;               % reflective rate of wall
R_head = 0.6;               % reflective rate of head
Z_wall = rho * c *((1+R_wall) / (1-R_wall));        % DECIDED BY absorbtion = 0.5
Z_head = rho * c *((1+R_head) / (1-R_head));

%% PML Layers
W = 100; 



% for n = 1 : 1 : W
%    Rx(n, :) = Rmax * ((W+1 - n) / W)^m;
%    Rx(nx+2*W+1-n, :) = Rx(n, :);
% end
% 
% for n = 1 : 1 : W
%    Ry(:, n) = Rmax * ((W+1 - n) / W)^m;
%    Ry(:, ny+2*W + 1 - n) = Ry(:, n);
% end


% 点音源の各定数
sig_freq = 1.0e3;           % 印可波形の周波数 [Hz]
sig_amp = 300.0;           % 印可波形の振幅
sig_x = (nx + 2 * W) / 2;                % 点音源の x 座標 [pixels] ( 0 < sig_x <= nx)
sig_y = (ny + 2 * W) / 2;                % 点音源の y 座標 [pixels] ( 0 < sig_y <= ny)


% initiate sound pressure
Vx = zeros(nx+2*W+1, ny+2*W);      % x方向粒子速度 [m/s]
Vy = zeros(nx+2*W,  ny+2*W+1);      % y方向粒子速度 [m/s]
P  = zeros(nx+2*W,  ny+2*W);      % 音圧 [Pa]
Px = zeros(nx+2*W,  ny+2*W); 
Py = zeros(nx+2*W,  ny+2*W); 
Precord = zeros(nx,1);


color_step = 50;
mymap = zeros(color_step * 2 + 1, 3);

for n = 1 : 1 : color_step
    mymap(n, :) = [n / color_step - 1 / color_step, n / color_step - 1 / color_step, 1];   
end
mymap(color_step + 1, :) = [1, 1, 1];

for n = 1 : 1 : color_step
    mymap(n + color_step  + 1, :) = [1, (color_step - n) / color_step, (color_step - n) / color_step];
end

%% load TSP

y = load("arbitrary_TSP.mat");
y = struct2cell(y);
y = 3000 * real(cell2mat(y));


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

%% 障害物のindexを取得
center_x = 200;
center_y = 100;
radius = 10;
field = zeros(nx, ny);
field2 = zeros(nx, ny);
field3 = zeros(nx, ny);

for xx = 1 : 1 : nx
    for yy = 1 : 1 : ny
       distance = sqrt((xx - center_x)^2 +  (yy - center_y)^2);
       if(distance > radius)
           field(xx, yy) = 0;
       elseif(distance < radius)
           field(xx, yy) = 1;
       end
    end
end

object_index = find(field == 1);

for xx = 2 : 1 : nx-1
   for yy = 2 : 1 : ny-1
      if (field(xx, yy) == 1)
         if (field(xx - 1, yy) == 0 && field(xx + 1, yy) == 1 && field(xx, yy - 1) == 1 && field(xx, yy + 1) == 0)
             field2(xx, yy) = 101;
         elseif (field(xx - 1, yy) == 1 && field(xx + 1, yy) == 1 && field(xx, yy - 1) == 1 && field(xx, yy + 1) == 0)
             field2(xx, yy) = 102;
         elseif (field(xx - 1, yy) == 1 && field(xx + 1, yy) == 0 && field(xx, yy - 1) == 1 && field(xx, yy + 1) == 0)
             field2(xx, yy) = 103;
         elseif (field(xx - 1, yy) == 0 && field(xx + 1, yy) == 1 && field(xx, yy - 1) == 1 && field(xx, yy + 1) == 1)
             field2(xx, yy) = 104;
         elseif (field(xx - 1, yy) == 1 && field(xx + 1, yy) == 0 && field(xx, yy - 1) == 1 && field(xx, yy + 1) == 1)
             field2(xx, yy) = 105;
         elseif (field(xx - 1, yy) == 0 && field(xx + 1, yy) == 1 && field(xx, yy - 1) == 0 && field(xx, yy + 1) == 1)
             field2(xx, yy) = 106;
         elseif (field(xx - 1, yy) == 1 && field(xx + 1, yy) == 1 && field(xx, yy - 1) == 0 && field(xx, yy + 1) == 1)
             field2(xx, yy) = 107;
         elseif (field(xx - 1, yy) == 1 && field(xx + 1, yy) == 0 && field(xx, yy - 1) == 0 && field(xx, yy + 1) == 1)
             field2(xx, yy) = 108;
         end
      end
   end
end

boundary101_index = find(field2 == 101);

boundary102_index = find(field2 == 102);

boundary103_index = find(field2 == 103);

boundary104_index = find(field2 == 104);

boundary105_index = find(field2 == 105);

boundary106_index = find(field2 == 106);

boundary107_index = find(field2 == 107);

boundary108_index = find(field2 == 108);

P(object_index) = 0;

%% start calculate sound field

for n = 1 : 1 : nmax
    
    Vx(2:end-1, :) = Vx(2:end-1, :) + (dt / rho) * ((P(1: end-1, :) - P(2:end, :)) / dx);
    Vy(:, 2:end-1) = Vy(:, 2:end-1) + (dt / rho) * ((P(:, 1: end-1) - P(:, 2:end)) / dx);
    
    
    %     update boundary
    Vx(1, :)   = - P(1, :)  / Z_wall;
    Vx(end, :) =   P(end, :) / Z_wall;
    Vy(:, 1)   = - P(:, 1)  / Z_wall;
    Vy(:, end) =   P(:, end) / Z_wall;
    
%     for index = 1 : 1 : length(boundary101_index)
%         [xx, yy] = ind2sub(size(P), boundary101_index(index));
%         Vx(xx, yy) = P(xx - 1, yy) / Z_head;
%         Vy(xx, yy + 1) = -P(xx, yy + 1) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary102_index)
%        [xx, yy] = ind2sub(size(P), boundary102_index(index));
%        Vy(xx, yy + 1) = -P(xx, yy + 1) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary103_index)
%         [xx, yy] = ind2sub(size(P), boundary103_index(index));
%         Vx(xx + 1, yy) = -P(xx + 1, yy) / Z_head;
%         Vy(xx, yy + 1) = -P(xx, yy + 1) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary104_index)
%         [xx, yy] = ind2sub(size(P), boundary104_index(index));
%         Vx(xx, yy) = P(xx - 1, yy) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary105_index)
%         [xx, yy] = ind2sub(size(P), boundary105_index(index));
%         Vx(xx + 1, yy) = -P(xx + 1, yy) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary106_index)
%         [xx, yy] = ind2sub(size(P), boundary106_index(index));
%         Vx(xx, yy) = P(xx - 1, yy) / Z_head;
%         Vy(xx, yy) = P(xx, yy - 1) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary107_index)
%         [xx, yy] = ind2sub(size(P), boundary107_index(index));
%         Vy(xx, yy) = P(xx, yy - 1) / Z_head;
%     end
%     
%     for index = 1 : 1 : length(boundary108_index)
%         [xx, yy] = ind2sub(size(P), boundary108_index(index));
%         Vx(xx + 1, yy) = -P(xx + 1, yy) / Z_head;
%         Vy(xx, yy) = P(xx, yy - 1) / Z_head;
%     end
    
    Px(:, :) = Px(:, :) + dt * ((kappa / dx) * (Vx(1:end-1, :) - Vx(2:end, :)));
    Py(:, :) = Py(:, :) + dt * ((kappa / dx) * (Vy(:, 1:end-1) - Vy(:, 2:end)));
    
    
    P(:, :) = (Px(:, :) + Py(:, :));

    
    % square cosine function
%     if n < (1.0/sig_freq)/dt,
%         sig(n + 1) = sig_amp * (1.0-cos(2.0*pi*sig_freq*n*dt))/2.0 * sin(2.0*pi*sig_freq*n*dt);
%         %         sig = sig_amp * sin(2.0*pi*sig_freq*n*dt);
%         P(sig_x,sig_y) = P(sig_x,sig_y) + sig(n + 1);
%     end

    % TSP
    if(n <= length(y))
        P(sig_x, sig_y) = P(sig_x, sig_y) + y(n);
    end
  
    
    
    
    
    
    
    Precord(n) = P(sig_x + 70, sig_y);
    
    
    
    
    
    
    
    
%     if(n == 2200)
%         figure(10);
%         clims = ([-5 5]);
%         imagesc((P(101:400, 101:300))', clims);     
%         axis equal; axis tight;             % 軸の調整
%         title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
%         xlabel('x'); ylabel('y');
%         colormap(mymap);
%         colorbar;
% %         saveas(gcf, './result/TSP_logweighted_snapshot.png');
%     end
    
    
    
    
    
%     mesh(abs(P)');
%     zlim([-5 3]);


    clims = [-5 5];
    imagesc((P)', clims);     % 行列の画像表示（正負の値を表示する場合）：↑いずれかのみ

    axis equal; axis tight;             % 軸の調整
	title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
    xlabel('x'); ylabel('y');
    colormap(mymap);
    colorbar;

    drawnow;  
%     
end

figure(1);
subplot(2,1,1);
plot(Precord);
xlabel('step'); 

TSP_freq = fft(Precord(801 : 4000));
subplot(2,1,2);
semilogx(1 : fs / length(TSP_freq) : fs/2, 20 * log10(abs(TSP_freq(1 : length(TSP_freq) / 2))));
xlabel('frequency[Hz]');

% saveas(gcf, './result/TSP_logweighted.png');


iTSP = load("arbitrary_TSP_inverse.mat");
iTSP = struct2cell(iTSP);
iTSP = 5000 * real(cell2mat(iTSP));

PULSE = conv(Precord, iTSP);
PULSE_freq = fft(PULSE);
figure(3);
subplot(2,1,1);
plot(PULSE);
subplot(2,1,2);
semilogx(1 : fs / length(PULSE_freq) : fs/2, 20 * log10(abs(PULSE_freq(1 : ceil(length(PULSE_freq) / 2) ) ) ) );









