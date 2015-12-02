% pos.m
% Positive function for the VITE model
% Antonio Ulloa
% Cognitive & Neural Systems
% Sun Oct  8 15:26:44 EDT 2000

function pos = pos(V)

     pos = (V>0).*V;
