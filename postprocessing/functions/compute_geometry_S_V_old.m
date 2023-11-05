clear all, close all, clc;

w1 = 30; % width
l1 = 10; % length
h1 = 10; % height

w2 = 26;
l2 = 10;
h2 = 10;

w3 = 20;
l3 = 10;
h3 = 10;

V1 = w1*l1*h1;
V2 = w2*l2*h2;
V3 = w3*l3*h3;

Vtot = V1+V2+V3

% 18 absorbing surfaces
S = 0;

% floor
S = S + w1*l1 + w2*l2 + w3*l3;

% ceiling
S = 2*S;

% walls
S = S + 2*(l1*h1 + l2*h2 + l3*h3);
S = S + w1*h1 + w3*h3;
S = S + (w1 - w2)*h2;
S = S + (w2 - w3)*h2;

S






