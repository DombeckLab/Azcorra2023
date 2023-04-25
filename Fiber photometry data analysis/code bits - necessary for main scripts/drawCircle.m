function [xunit, yunit] = drawCircle(x,y,r)
% draw a circle with center at coordinates (x,y) and radius r
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end