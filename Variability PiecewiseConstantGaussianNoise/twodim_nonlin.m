delT = 1;
x = [ 0; -2; 0 ];
u = [ -1; 0.1 ];

XB = [5, 5];

Q = [ 0.1 0    0
      0    0.1 0
      0    0     0.0002];
R = [ 0.01 0
      0  0.00001 ];
   
%% Correct start is [0; -2; 0].  Also try [0; 0; 0]    
xhat = [0; -2; 0 ];
P = 10*Q;

for i=1:100
   % generate true state
   w = gennormal([0;0;0], Q);          % generate process noise
   x = robot(x, u); %+ w;          % update state
   x(1:2)
   %plot(x(1),x(2),'m+');
   plotrobot(x,'b');
   hold on;
   plot(XB(1),XB(2),'ro');
   axis([-10 10 -10 10]);
   pause;   

   % generate measurement
   v = gennormal([0;0], R);
   z = beacon(x,XB) + v;
  
   % predict
   [xpred, Ppred] = predict_nonlin('robot', xhat, u, P, Q);
   plot(xpred(1), xpred(2),'+g');
   plot_covar(xpred(1:2),Ppred(1:2,1:2),'g');
   pause;

   % update
   [xhat, P] = update_nonlin(xpred, Ppred, z, 'beacon', XB, R);
   plot(xhat(1),xhat(2),'+r');
   plot_covar(xhat(1:2),P(1:2,1:2),'r');
   pause;

   %hold off;   
end
