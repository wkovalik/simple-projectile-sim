clear; clc; close all;


X1 = zeros(2, 10000);
X2 = zeros(2, 10000);
for i = 1:10000
    X1(:, i) = randn(2, 1);
    X2(:, i) = [randn(); randn()];
end

plot(X1(1, :), X1(2, :), '.')
hold on
plot(X2(1, :), X2(2, :), '.')
hold off