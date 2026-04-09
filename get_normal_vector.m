function normal = get_normal_vector(plane_2D)
    % 返回指向体积外侧的单位法向量（outward normal）
    switch plane_2D
        case "yz_x=-2N"
            normal = [-1, 0, 0];  % x = -2N 面，外法线指 -x
        case "yz_x=2N"
            normal = [1, 0, 0];   % x = +2N 面，外法线指 +x
        case "xz_y=-2N"
            normal = [0, -1, 0];  % y = -2N 面，外法线指 -y
        case "xz_y=2N"
            normal = [0, 1, 0];   % y = +2N 面，外法线指 +y
        case "xy_z=-2N"
            normal = [0, 0, -1];  % z = -2N 面，外法线指 -z
        case "xy_z=2N"
            normal = [0, 0, 1];   % z = +2N 面，外法线指 +z
        otherwise
            error('Unknown plane: %s', plane_2D);
    end
end
