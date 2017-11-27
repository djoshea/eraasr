function [x, idxFirst, idxLast] = infillNanEdgesWithLastSample(x)
% operates along dimension 1, takes the first non-nan sample and fills in
% backwards with that value, then the last non-nan sample and fills in
% forwards with that value

szX = size(x);
x = x(:, :);

[idxFirst, idxLast] = deal(nan([1 szX(2:end)]));

for c = 1:size(x, 2)
    mask = ~isnan(x(:, c));
    if ~any(mask)
        continue;
    end
    first = find(mask, 1, 'first');
    last = find(mask, 1, 'last');
    x(1:first-1, c) = x(first, c);
    x(last+1:end, c) = x(last, c);
    
    idxFirst(c) = first;
    idxLast(c) = last;
end

x = reshape(x, szX);

    
    

