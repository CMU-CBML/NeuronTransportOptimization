function [ P ] = RotateAroundAxis( Q, axis, angle )
%RotateAroundAxis: rotate the geometry around defined axis with angle
axis=axis/norm(axis);

cos_a = cos(angle);
sin_a=sin(angle);

R= [cos_a+axis(1)^2*(1-cos_a), axis(1)*axis(2)*(1-cos_a)-axis(3)*sin_a, axis(1)*axis(3)*(1-cos_a)+axis(2)*sin_a;
    axis(2)*axis(1)*(1-cos_a)+axis(3)*sin_a,cos_a+axis(2)^2*(1-cos_a), axis(2)*axis(3)*(1-cos_a)-axis(1)*sin_a;
    axis(3)*axis(1)*(1-cos_a)-axis(2)*sin_a, axis(3)*axis(2)*(1-cos_a)+axis(1)*sin_a,cos_a+axis(3)^2*(1-cos_a)];

P=Q*R';
end

