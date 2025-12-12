clear all, close all, clc;

out_filename = 'out_0.txt';
room_filename = "hall.txt";
dh = 0.5;

A = load(out_filename);
N = size(A,1);

[Lx_int, Ly_int, Lz_int] = calculate_grid_size(room_filename,dh);
Lx = double(Lx_int);
Ly = double(Ly_int);
Lz = double(Lz_int);

% Reshape the data into a 3D array
pressure_values = zeros(N,Lx,Ly,Lz);

for iT = 1:N
    pressure_values(iT,:,:,:) = reshape(A(iT,:), Lx, Ly, Lz);
end

%%
close all;

for iT = 2:1:5
    figure()
    p = squeeze(pressure_values(iT,:,:,:));
    
    slice(p, Ly*5/6, Lx/2, Lz*3/10);  % Create a slice plot
    
    shading interp
    % title(['Pressure Distribution at Time Instant ', num2str(time_instants(iT))]);
    daspect([1, 1, 1]); % Equal aspect ratio for x, y, and z dimensions
    xlabel('y')
    ylabel('x')
    zlabel('z')
    clim([-1 1])
    colormap jet
    % Change font size
    set(gca, 'FontSize', 24)
end
drawnow;