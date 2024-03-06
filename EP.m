classdef EP < handle
    properties
        constellation {mustBeNumeric}
        constellation_len {mustBeNumeric}
        min_var {mustBeNumeric} = eps       % the default minimal variance is 2.2204e-16
        iter_num = 10                       % maximal iteration
        iter_diff_min = eps;                % the minimal difference between 2 adjacent iterations
    end
end