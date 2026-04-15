clear; clc; close all;

tic
x = rand(6, 1);
for i = 1:1000000
    new_x = rand(6, 1);
    x = new_x;
end
toc

tic
x = rand(6, 1);
for i = 1:1000000
    new_x = rand(6, 1);
    x(1:6) = new_x(1:6);
end
toc

tic
x = rand(6, 1);
for i = 1:1000000
    new_x = rand(6, 1);
    x(1) = new_x(1);
    x(2) = new_x(1);
    x(3) = new_x(1);
    x(4) = new_x(1);
    x(5) = new_x(1);
    x(6) = new_x(1);
end
toc
