function resolved = resolve_contrasts(contrasts, conditionNames)
% Expand named interaction terms into the source model's condition order.
resolved = struct('name', {}, 'weights', {});
for k = 1:numel(contrasts)
    c = contrasts(k);
    if isfield(c, 'conditions') && ~isempty(c.conditions)
        terms = string(c.conditions);
        weights = c.weights;
        assert(isnumeric(weights) && isreal(weights) && isvector(weights) && ...
            numel(weights) == numel(terms) && all(isfinite(weights)) && any(weights), ...
            'Contrast "%s" needs one finite weight per named condition and a nonzero effect.', c.name);
    else
        assert(isfield(c, 'positive') && isfield(c, 'negative') && ...
            isscalar(string(c.positive)) && isscalar(string(c.negative)), ...
            'Contrast "%s" needs positive and negative condition names, or conditions and weights.', c.name);
        terms = [string(c.positive), string(c.negative)];
        weights = [1 -1];
    end
    assert(numel(unique(terms)) == numel(terms), 'Repeated condition in contrast "%s".', c.name);
    ordered = zeros(1, numel(conditionNames));
    for j = 1:numel(terms)
        match = conditionNames == terms(j);
        assert(nnz(match) == 1, 'Contrast "%s": expected one match for "%s". Available conditions: %s.', ...
            c.name, terms(j), strjoin(conditionNames, ', '));
        ordered(match) = weights(j);
    end
    resolved(k).name = char(c.name);
    resolved(k).weights = ordered;
end
end
