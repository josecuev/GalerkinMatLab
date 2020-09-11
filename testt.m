x = rand(5,1);
y = rand(5,1);
[r,t] = cart2pol(x,y);
fprintf('\n\n%11s%11s%11s%11s\n','x', 'y', 'r', 'theta');
fprintf(' %10.2f %10.2f %10.2f %10.2f\n',[x,y,r,t].');
fprintf('\n\n')