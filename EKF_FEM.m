function[x_final,P,zhat]=EKF_FEM(x_in,Z,Options,P,R)

XMap = Options.XMap;
YMap = Options.YMap;
Input_geometry = Options.inputPara;

%x_in=[E1,E2,E3]
%inial condition
xi=length(x_in);
xo=x_in;
Phi=eye(xi);
I3=eye(xi);
R=R*eye(length(Z));

qd=0.000002; %Discrete process noise covariacne 
Q=qd*eye(xi);


x=xo;
% calculating the deflection corresponding to the guessed E
zhat=h_z(x,Input_geometry,XMap,YMap);

% calculating the Jacobian (gradient) at the current guessed E
[H]=Jacobian_H(x,Input_geometry,XMap,YMap);

%computing a priori covariacne matrix      
P=Phi*P*Phi'+Q;
%computing the Kalman Gain
K=P*H'/[H*P*H'+R];
%2.measurment update
%conditioning the predicted estimate on the measurment
x=x(:)+K*(Z-zhat);
%computing the a posteriori covaricance matrix
P=(I3-K*H)*P;
x_final=x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Observility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OB=obsv(Phi,H);
%Observ=rank(OB);

%plot(reference_point,U1_out)
%grid on
%title('Displacement')
%xlabel('length (in)')
%ylabel('displacement (in) to x direction')
%set(gca,'PlotBoxAspectRatio',[1 1 1])



function [dhdx]=Jacobian_H(x,Input_geometry,XMap,YMap)

E = x;
dE = E*0.01;

% Make sure the perturbation is not zero
dE(dE<0.01) = 0.01;

E0 = repmat(E,[1,length(E)]);
E1 = E0-diag(dE);
E2 = E0+diag(dE);

E1E2=[E1,E2];
H = h_z(E1E2,Input_geometry,XMap,YMap);
H1 = H(:,1:length(E));
H2 = H(:,length(E)+1:end);

dx = ones(size(H1,1),1)*dE';
dhdx = (H2-H1)/2./dx;