function plotReachTrajectories(X,Y,Z,X_from_under,reachTrajTimes,smoobin)

if ~isempty(smoobin)
    X=smooth(X,smoobin);
    Y=smooth(Y,smoobin);
    Z=smooth(Z,smoobin);
    X_from_under=smooth(X_from_under,smoobin);
end

figure(); 
plot3(X,Y,Z); 
hold all; 
plot3(X(1:50),Y(1:50),Z(1:50)); 
xlabel('X'); ylabel('Y'); zlabel('Z');