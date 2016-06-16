function points = lin_tri(m1,m2,P1,P2,N)
%=======================================================================
% The function implements linear triangulation method
% which is a direct analog of DLT as described by Richard.H and Andrew.Z 
% in 'Multiple View Geometry' in section 12.2, page 312.
% input --> coordinates of corners in the image pair and N = no. of matches
% output --> 3d Reconstruction points
%=======================================================================
X = ones(4, N);
%Combinging the measurement x = PX. x' = P'X, from the image to form AX = 0
for i = 1:N
    A =[m1(1,i)*P1(3,1:3) - P1(1,1:3);...
        m1(2,i)*P1(3,1:3) - P1(2,1:3);...
        m2(1,i)*P2(3,1:3) - P2(1,1:3);...
        m2(2,i)*P2(3,1:3) - P2(2,1:3)];
    
    b = -[m1(1,i)*P1(3,4) - P1(1,4);...
        m1(2,i)*P1(3,4) - P1(2,4);...
        m2(1,i)*P2(3,4) - P2(1,4);...
        m2(2,i)*P2(3,4) - P2(2,4)];
    %the unit singular vector corresponding to the smallest singular value of A
    [u,s,v] = svd(A);
    B = u'*b;
    y = B(1:3) ./ (diag(s));
    X(1:3,i) = v*y;
end
X = X./repmat(X(4,:),4,1);
X = X(1:3, :);
points = X';
end