function plotReachTrajectories(X,Y,Z,X_from_under,reachTrajTimes,smoobin)

if ~isempty(smoobin)
    for i=1:size(X,1)
        X(i,:)=smooth(X(i,:)',smoobin);
        Y(i,:)=smooth(Y(i,:)',smoobin);
        Z(i,:)=smooth(Z(i,:)',smoobin);
        X_from_under(i,:)=smooth(X_from_under(i,:)',smoobin);
    end
end

figure(); 
for i=1:size(X,1)
    plot3(X(i,:),Y(i,:),Z(i,:),'Color','k'); 
    hold all; 
    plot3(X(i,1:50),Y(i,1:50),Z(i,1:50),'Color','b'); hold on; 
end
xlabel('X'); ylabel('Y'); zlabel('Z');