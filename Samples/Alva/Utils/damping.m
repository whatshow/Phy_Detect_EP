function [out] = damping(old,new,m)
% compute the dumping of b and c by a factor m
out = m .* old + (1 - m) .* new;

end

