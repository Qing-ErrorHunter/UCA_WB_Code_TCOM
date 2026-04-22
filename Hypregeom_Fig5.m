clear all;
clc;
close all;

%% hypergeometric function 1F2
a = 1/2;
b = [1,3/2];

C = linspecer(9);
z_max = 20;
z_list = 0:0.01:z_max;
g = hypergeom(a,b,-z_list.^2/4);
figure;
hold on;
plot(z_list, g,'-','LineWidth',1.5,'color', C(1,:));
xlabel('x');
ylabel('Value');
grid on;
box on;

%% hypergeometric function 2F3
a2 = [1/2,1/2];
b2 = [1,3/2,3/2];

g2 = hypergeom(a2,b2,-z_list.^2/4);
plot(z_list, g2,'-','LineWidth',1.5,'color', C(2,:));
legend('$_1F_2(1/2;1,3/2;-x^2/4)$','$_2F_3(1/2,1/2;1,3/2,3/2;-x^2/4)$','Interpreter','latex');
