function [F,residual_mean] = fit_fundamental(m)
%============================================================================
% Function computes the fundamental matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used and is based on Peter Kovesi's code
%(http://www.peterkovesi.com/matlabfns/)
% Copyright (c) 2002-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% Code used by Syed Ahmed for Assignment 4 of CAP 6419 at UCF
%============================================================================
    F = zeros(3,3); % initializing the final output...
    threshold = 2;
    nom = size(m, 1); % number of matches...
	m1 = m(:,1:2); % spliting the set of points...
	m2 = m(:,3:4); % ... based on image 1 and 2...
	% Normalising the first set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
	u1 = mean(m1(:, 1));
    m1(:,1) = m1(:,1) - u1;
    u2 = mean(m1(:, 2));
    m1(:,2) = m1(:,2) - u2;
    dist = sqrt(m1(:,1).^2 + m1(:,2).^2);
    dist_mean = mean(dist);
    scale = 2/dist_mean;
    m1 = m1*scale;
    T1 = [scale    0     -scale*u1
           0    scale   -scale*u2
           0      0         1   ];
	m(:,1:2) = m1;
	% Normalising the second set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2).
	u1 = mean(m2(:, 1));
    m2(:,1) = m2(:,1) - u1;
    u2 = mean(m2(:, 2));
    m2(:,2) = m2(:,2) - u2;
    dist = sqrt(m2(:,1).^2 + m2(:,2).^2);
    dist_mean = mean(dist);
    scale = 2/dist_mean;
    m2 = m2*scale;
    T2 = [scale    0     -scale*u1
           0    scale   -scale*u2
           0      0         1   ];
	m(:,3:4) = m2;
	
    in_pts = [1:nom]'; %inlier points....
    % buitding the constraint matrix...
	A = [];
	for i = 1:size(in_pts,1)
        u = m(in_pts(i), 1);
        v = m(in_pts(i), 2);
        u1 = m(in_pts(i), 3);
        v1 = m(in_pts(i), 4); 
        A = [A; u1*u u1*v u1 v1*u v1*v v1 u v 1];
    end
    [U,S,V]=svd(A); 
	% Extract fundamental matrix from the column of V corresponding to
    % smallest singular value.
    F = reshape(V(:,end), 3, 3)';
	% Enforce constraint that fundamental matrix has rank 2 by performing
    % a svd and then reconstructing with the two largest singular values.
    [U,S,V] = svd(F,0);
    F = U*diag([S(1,1) S(2,2) 0])*V';
	% Denormalising...
    F = T2'* F * T1;
	residual = [];
    for i = 1:nom
        temp = F * [m(i, 1) m(i, 2) 1]';
        residual = [residual; abs([m(i, 3),m(i, 4), 1]*temp)];
    end
	residual;
	residual_mean = mean(residual);
end