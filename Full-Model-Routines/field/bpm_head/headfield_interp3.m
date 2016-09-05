function h_interp = headfield_interp3(x, y, z, x_data, y_data, z_data, h_data)
h_interp = interp3(x_data, y_data, z_data, h_data, x, y, z, 'cubic'); %'spline'); % gets out of memory message with spline
end
