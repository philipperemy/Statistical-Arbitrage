function s = strstartswith(a, b)

s = (numel(a)>=numel(b)) && all(a(1:min(numel(a),numel(b)))==b(1:min(numel(a),numel(b))));