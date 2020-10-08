  
function r = spaceTimeCovarianceExpGeom(di1,dj1,t1,di2,dj2,t2,thetas,thetad,thetat)
dist = sqrt(abs(di1-di2)^2+abs(dj1-dj2)^2);
distt = abs(t1-t2);

distSpaceTime = sqrt((dist/thetad)^2 + (distt/thetat)^2);

r = thetas * exp(-distSpaceTime);
