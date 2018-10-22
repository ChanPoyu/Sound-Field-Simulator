%% FDTD PML

clear all;
close all;


% ���f���̊e�萔
nx = 300;                   % X �����̋�ԃZ���� [pixels]
ny = 200;                   % Y �����̋�ԃZ���� [pixels]
dx = 10.0e-3;               % ��ԍ��� [m]
fs = 48000;
dt = 1/fs;                  % ���ԍ��� [s]
nmax = 10000;                % �v�Z�X�e�b�v�� [��]
savestep = 10;              % �ۑ��Ԋu�i�X�e�b�v���j
f = linspace(0 ,  24000, nmax / 2);


% �}���̊e�萔
rho = 1.293;                % �}���̖��x�� [kg/m^3]
kappa = 142.0e3;            % �}���̑̐ϒe������ [Pa]
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


% �_�����̊e�萔
sig_freq = 1e3;           % ��g�`�̎��g�� [Hz]
sig_amp = 100.0;           % ��g�`�̐U��
sig_x = (nx + 2 * W) / 2;                % �_������ x ���W [pixels] ( 0 < sig_x <= nx)
sig_y = (ny + 2 * W) / 2;                % �_������ y ���W [pixels] ( 0 < sig_y <= ny)


% initiate sound pressure
Vx = zeros(nx+2*W+1, ny+2*W);      % x�������q���x [m/s]
Vy = zeros(nx+2*W,  ny+2*W+1);      % y�������q���x [m/s]
P  = zeros(nx+2*W,  ny+2*W);      % ���� [Pa]
Px = zeros(nx+2*W,  ny+2*W); 
Py = zeros(nx+2*W,  ny+2*W);
Psnapshot = zeros(nx+2*W,  ny+2*W);
Precord = zeros(1,nmax);

%% mic position
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
y = 30 * real(cell2mat(y));


%% Gaussian pulse

% a = 3;
% for xx = 1 : 1 : nx+2*W
%     for yy = 1 : 1 : ny+2*W
%         r_sqr = (xx - 150)^2 + (yy - 200)^2;
%         if r_sqr < 10000
%             P(xx, yy) = sig_amp * exp(- r_sqr / a^2);
%         end
%     end
% end

% P(150 , 200) = 300;

F(10000) = struct('cdata',[],'colormap',[]);
writerObj = VideoWriter('Primary.avi');
writerObj.FrameRate = 10;

%% start calculate 

for n = 1 : 1 : nmax
    
    Vx(2:end-1, :) = Vx(2:end-1, :) - (dt / rho) * Vx(2:end-1, :) .* Rx(2:end, :) - (dt / rho) * ((P(1: end-1, :) - P(2:end, :)) / dx);
        Vy(:, 2:end-1) = Vy(:, 2:end-1) - (dt / rho) * Vy(:, 2:end-1) .* Ry(:, 2:end) - (dt / rho) * ((P(:, 1: end-1) - P(:, 2:end)) / dx);
        
        
        Px(:, :) = Px(:, :) - dt * times(Px(:, :), Rx(:, :)) / rho - dt * ((kappa / dx) * (Vx(1:end-1, :) - Vx(2:end, :)));
        Py(:, :) = Py(:, :) - dt * times(Py(:, :), Ry(:, :)) / rho - dt * ((kappa / dx) * (Vy(:, 1:end-1) - Vy(:, 2:end)));
        
        
        P(:, :) = (Px(:, :) + Py(:, :));

    
    
    
    
    
    
    % square cosine function
%         sig = sig_amp * sin(2.0*pi*sig_freq*n*dt);
% 
%         P(150,200) = P(150,200) + sig;

%     TSP
    if(n <= length(y))
%         P(150, 200) = P(150, 200) + rho * c ^ 2 * dt * y(n) /  (1e-10 + sqrt(Vx(150,200)^2 + Vy(150,200)^2)) ;
%         P(150, 200) = P(150, 200) + rho * c ^ 2 * dt * y(n);
        P(150, 200) = P(150, 200) +  + y(n);
    end
  
    
    
    
    
    
    
    Precord(n) = P(350, 200);
    
    
    
    
    
    
    
%     
%     if(n == 2200)
%         figure(10);
%         clims = ([-5 5]);
%         imagesc((P(101:400, 101:300))', clims);     
%         axis equal; axis tight;             % ���̒���
%         title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
%         xlabel('x'); ylabel('y');
%         colormap(mymap);
%         colorbar;
%         plotmic = plot(mic_position(:, 1) - 100, mic_position(:, 2) - 100, 'k');
%         plotmic.Marker = '*';
%         plotmic.LineWidth = 2;
% %         saveas(gcf, './result/TSP_logweighted_snapshot.png');
%     end
    
    
%     if (n > 1000)
    
%         clims = [-5 5];
%         imagesc((P(101 : 400, 101 : 300))', clims);     % �s��̉摜�\���i�����̒l��\������ꍇ�j�F�������ꂩ�̂�
% 
%         axis equal; axis tight;             % ���̒���
%         title(sprintf('%7.3f ms  [ %4d/%4d ]   ([Ctrl]-[c] to STOP)', n*dt*1.0e3, n, nmax));
%         xlabel('x'); ylabel('y');
%         colormap(mymap);
%         colorbar;
%         hold on;
%         plotmic = plot(mic_position(:, 1) - 100, mic_position(:, 2) - 100, 'k');
%         plotmic.Marker = '*';
%         plotmic.LineWidth = 2;
% 
%         F(n) = getframe(gcf);
        drawnow;
    if(rem(n, 80) == 0)
        close all;
    end
    
    if(n == 50)
        Psnapshot(:, :) = P(:, :);
    end
    
end

energy = sum(sum(Psnapshot));
 
% open(writerObj);
% for n = 1 : 1 : length(F)
%    frame = F(n);
%    writeVideo(writerObj, frame);
% end
% close(writerObj);

Precord = Precord / max(Precord);

figure(1);
subplot(1,2,1);
plot(Precord);
title('[FDTD]Sine Square in frequency domain');
xlabel('step');
axis tight;


Y = fft(Precord);
subplot(1,2,2);
semilogx(f, 10 * log10(abs(Y(1 : (end)/2))));
title('[FDTD]Sine Square in frequency domain');
xlabel('[Hz]');

axis tight;


