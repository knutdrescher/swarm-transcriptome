function pos = getPos(name)
    match = regexp(name, ...
        '(?<=^x)(?<pos_x>-?\d+)y(?<pos_y>-?\d+)(?=.*$)','names');

    x = str2double(match.pos_x);
    y = str2double(match.pos_y);
    pos = sqrt(x^2 + y^2);
end