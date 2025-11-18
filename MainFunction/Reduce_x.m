function [tmp] = Reduce_x(tmp, t, inds_t, x, cellInp, prm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:prm.m

    % 1. Update Pt(t)
    tmp{i}.Pt(t) = tmp{i}.Pt(t) - tmp{i}.Px;

    % 2. Pre-calculate Pyx(:,x)
    Pyx_x = cellInp{i}.Pyx(:, x);

    % 3. Remove x from inds_t (efficiently)
    inds_t(inds_t == x) = [];

    % 4. Incremental update of Py_t(:,t)
    %   - Subtract the contribution of x from the *total* sum represented by
    %     Py_t(:,t) * Pt(t) *before* removing x.
    %   - Then, divide by the *new* Pt(t).

    % Check for division by zero *before* calculations
    if abs(tmp{i}.Pt(t)) < 1e-12
        tmp{i}.Py_t(:, t) = zeros(size(cellInp{i}.Pyx, 1), 1); % Or handle differently
    else
        % Calculate the updated Py_t(:,t)
        tmp{i}.Py_t(:, t) = (tmp{i}.Py_t(:, t) * (tmp{i}.Pt(t) + tmp{i}.Px) - Pyx_x) / tmp{i}.Pt(t);

        % Handle small values in Py_t(:,t)
        f_small_values = abs(tmp{i}.Py_t(:, t)) < 1e-12;
        tmp{i}.Py_t(f_small_values, t) = 0;
    end
end

end

