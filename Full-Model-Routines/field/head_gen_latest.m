function [hx_data hy_data hz_data x_vec y_vec z_vec] = head_gen_latest(headfield, z_extrap)

numx = 601; % number of x points: cross track
numy = 601; % number of y points: down track

i_vec = 2:7;
len_i_vec = length(i_vec);

% x_data = zeros(numx,numy,len_i_vec);
% y_data = zeros(numx,numy,len_i_vec);
z_data = zeros(numx,numy,len_i_vec);

hx_data = zeros(numx,numy,len_i_vec);
hy_data = zeros(numx,numy,len_i_vec);
hz_data = zeros(numx,numy,len_i_vec);
% h_eff = zeros(numx,numy,len_i_vec);

% y-coordinate is down track, x-coordinate is cross track
x_data = reshape(headfield(1:numx*numy, 1), numx, numy)';
y_data = reshape(headfield(1:numx*numy, 2), numx, numy)';
for i=1:len_i_vec
    %     x_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 1), numx, numy)';
    %     y_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 2), numx, numy)';
    z_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 3), numx, numy)';
    hx_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 4), numx, numy)';
    hy_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 5), numx, numy)';
    hz_data(:,:,i) = reshape(headfield(1 + (i-1)*numx*numy:i*numx*numy, 6), numx, numy)';
    %     h_eff(:,:,i) = headfield_eff(hx_data(:,:,i), hy_data(:,:,i), hz_data(:,:,i));
end

x_vec = x_data(1,:)'; % column vector 0:600
y_vec = y_data(:,1); % column vector 0:600
z_vec = reshape(z_data(1:1,1:1,:),1,length(z_data(1:1,1:1,:)))';

if (z_extrap ~= -1)
    
    disp('carrying out extrapolation');
    hx_extrap = extrap_h_data(z_extrap, x_vec, y_vec, z_vec, hx_data);
    hy_extrap = extrap_h_data(z_extrap, x_vec, y_vec, z_vec, hy_data);
    hz_extrap = extrap_h_data(z_extrap, x_vec, y_vec, z_vec, hz_data);

    nhx_data = zeros(numx,numy,len_i_vec+1);
    nhy_data = zeros(numx,numy,len_i_vec+1);
    nhz_data = zeros(numx,numy,len_i_vec+1);
    nz_vec = zeros(1,length(z_vec)+1);

    nz_vec(1) = z_extrap;
    nhx_data(:,:,1) = hx_extrap;
    nhy_data(:,:,1) = hy_extrap;
    nhz_data(:,:,1) = hz_extrap;
    for i=2:len_i_vec+1
        nz_vec(i) = z_vec(i-1);
        nhx_data(:,:,i) = hx_data(:,:,i-1);
        nhy_data(:,:,i) = hy_data(:,:,i-1);
        nhz_data(:,:,i) = hz_data(:,:,i-1);
    end

    z_vec = nz_vec;
    hx_data = nhx_data;
    hy_data = nhy_data;
    hz_data = nhz_data;
else
    disp('no extrapolation');
    % do nothing
end

end