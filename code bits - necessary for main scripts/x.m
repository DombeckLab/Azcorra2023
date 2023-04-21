function xaxis = x(data, varargin)
try
	xaxis = 1/data.hz * (1:length(data.vel));
catch
    try
        xaxis = 1/varargin{1} * (1:length(data));
    catch
        warning('hz = 100')
        xaxis = 1/100 * (1:length(data));
    end
end