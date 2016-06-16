%=====================================================================
%The following code is based on the sample code provided by Dr.Foroosh
%for the following tasks : (1) Displaying image pair side-by-side with
%matches; (2) finding epipolar lines and (3) display points and segments
%epipolar lines found in (2)
%=====================================================================

%=============================Step 1==================================
clear all;
close all;
while 1==1
	opt = input('what pair of images would you like to input?\nPress:\nH for house pair and\nL for library pair of images\n','s');
	if opt=='H' || opt=='h'
		I1 = imread('house1.jpg');
		I2 = imread('house2.jpg');
		P1 = load('house1_camera.txt');
		P2 = load('house2_camera.txt');
		matches = load('house_matches.txt');
		break;
	elseif opt=='L' || opt=='l'
		I1 = imread('library1.jpg');
		I2 = imread('library2.jpg');
		P1 = load('library1_camera.txt');
		P2 = load('library2_camera.txt');
		matches = load('library_matches.txt');
		break;
	else
		fprintf('Only H and L keys are allowed, please select the right choice\n')
	end
end
N = size(matches,1); % ttal features...
%=============================Step 2==================================
%% display two images side-by-side with matches
%% this code is to help you visualize the matches, 
%% you don't need to use it to produce the results for the assignment
%% display second image with epipolar lines reprojected 
%% from the first image
figure
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');
title('Matches For The Image Pair')
%pause;
hold off;
%=============================Step 3==================================
% first, fit fundamental matrix to the matches
[F,residual_mean] = fit_fundamental(matches); % this is the function that you need 
%% to build to call Peter Kovesi's 8-point algorithm implementation
fprintf('the mean the residual error is %f \n',residual_mean)
%----------------------------------------------------------
L = (F * [matches(:,1:2) ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% the mean squared distance in pixels between points in both images and
% the corresponding epipolar lines.
msd = zeros(N,1);
for i = 1:N
    msd(i,1) = sqrt(dist2(closest_pt(i,:), matches(i,1:2)));
end

fprintf('the mean squared distance is %f \n',mean(msd))

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
%clf;
figure
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');
title('Points and Segments of Epipolar Lines ')
hold off;

%=============================Step 4a==================================
% In this part of the code we implement the linear triangulation method
% which is a direct analog of DLT as described by Richard.H and Andrew.Z 
% in 'Multiple View Geometry' in section 12.2, page 312.
points1 = lin_tri(matches(:,1:2)',matches(:,3:4)',P1,P2,N);
%=============================Step 4b==================================
% In this part of the code we implement the optimal triangulation method
% which is a direct analog of DLT as described by Richard.H and Andrew.Z 
% in 'Multiple View Geometry' in section 12.5, page 318, algorithm 12.1
[x1_,x2_] = opt_tri(P1, P2, matches(:,1:2)', matches(:,3:4)', F,N);
points2 = lin_tri(x1_,x2_,P1,P2,N);
%=============================Step 4c==================================
% Finding the two camera centers to compare the reconstructed
% 3D points from the above two algorithms...
[u,s,v] = svd(P1);
center1 = v(:,end);
center1 = center1(1:3) ./ center1(4);
[u,s,v] = svd(P2);
center2 = v(:,end);
center2 = center2(1:3) ./ center2(4);
% to display the centers and reconstructed points...
% for linear triangulation method
figure
title('3D Reconstruction with Linear Triangulation Method')
plot3(points1(:,1),points1(:,2),points1(:,3),'ok');
hold on
plot3(center1(1),center1(2),center1(3),'Xb');
plot3(center2(1),center2(2),center2(3),'Xr');
legend('Reconstructed Points','Camera Center1','Camera Center2')
axis equal
%%-----------------------------------------------------------------
%% using Alan Jennings's code to create a video of the Reconstruction 
%% view rotating viewZ
%%-----------------------------------------------------------------
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true; 
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'Lin_Tri_2',OptionZ)
hold off
figure
title('3D Reconstruction with Optimal Triangulation Method')
plot3(points2(:,1),points2(:,2),points2(:,3),'ok');
hold on
plot3(center1(1),center1(2),center1(3),'Xb');
plot3(center2(1),center2(2),center2(3),'Xr');
axis equal
legend('Reconstructed Points','Camera Center1','Camera Center2')
%%-----------------------------------------------------------------
%% using Alan Jennings's code to create a video of the Reconstruction 
%% view rotating viewZ
%%-----------------------------------------------------------------
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'Opt_Tri_2',OptionZ)
hold off
%=============================Step 5==================================
% To Project the points on image and calcualte the the mean square errors
Lin1 = zeros(N, 3);
for i = 1 : N
    Lin1(i, :) = P1 * [points1(i, :), 1]';
    Lin2(i, :) = P2 * [points1(i, :), 1]';
end

for i = 1 : N
    Lin3(i, :) = Lin1(i, 1:2) ./ Lin1(i, 3);
    Lin4(i, :) = Lin2(i, 1:2) ./ Lin2(i, 3);
end
% for optimal triangulation method
Opt1 = zeros(N, 3);
for i = 1 : N
    Opt1(i, :) = P1 * [points2(i, :), 1]';
    Opt2(i, :) = P2 * [points2(i, :), 1]';
end

for i = 1 : N
    Opt3(i, :) = Opt1(i, 1:2) ./ Opt1(i, 3);
    Opt4(i, :) = Opt2(i, 1:2) ./ Opt2(i, 3);
end

figure
subplot(1,2,1)
imshow(I1)
hold on
fig1 = plot( matches(:, 1), matches(:, 2), 'xc' ) % given points...
fig2 = plot( Lin3(:, 1), Lin3(:, 2), 'xy' ) % reconstructed points...
hold off
%legend('Given Points','Reconstructed Points')
subplot(1,2,2)
imshow(I2)
hold on
fig3 = plot( matches(:, 3), matches(:, 4), 'xc' ), % given points...
fig4 = plot( Lin4(:, 1), Lin4(:, 2), 'xy' ), % reconstructed points...
legend('Given Points','Reconstructed Points')
hold off
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Comaprision of Given and Reconstructed Points on Image-pair for the Linear triangulation Method ','HorizontalAlignment','center','VerticalAlignment', 'top');
figure
subplot(1,2,1)
imshow(I1)
hold on
fig1 = plot( matches(:, 1), matches(:, 2), 'xc' ) % given points...
fig2 = plot( Opt3(:, 1), Opt3(:, 2), 'xy' ) % reconstructed points...
hold off
%legend('Given Points','Reconstructed Points')
subplot(1,2,2)
imshow(I2)
hold on
fig3 = plot( matches(:, 3), matches(:, 4), 'xc' ), % given points...
fig4 = plot( Opt4(:, 1), Opt4(:, 2), 'xy' ), % reconstructed points...
legend('Given Points','Reconstructed Points')
hold off
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Comaprision of Given and Reconstructed Points on Image-pair for the Optimal triangulation Method ','HorizontalAlignment','center','VerticalAlignment', 'top');
% MSE for Linear triangulation Method
mse_Lin = zeros(2, 1);
mse_Lin(1,:) = sum( ( matches(:, 1) - Lin3(:, 1) ) .^ 2 );
mse_Lin(1,:) = mse_Lin(1,:) + sum( ( matches(:, 2) - Lin3(:, 2) ) .^ 2 );
mse_Lin(2,:) = sum( ( matches(:, 3) - Lin4(:, 1) ) .^ 2 );
mse_Lin(2,:) = mse_Lin(2,:) + sum( ( matches(:, 4) - Lin4(:, 2) ) .^ 2 );
% MSE for Optimal triangulation Method
mse_Opt = zeros(2, 1);
mse_Opt(1,:) = sum( ( matches(:, 1) - Opt3(:, 1) ) .^ 2 );
mse_Opt(1,:) = mse_Opt(1,:) + sum( ( matches(:, 2) - Opt3(:, 2) ) .^ 2 );
mse_Opt(2,:) = sum( ( matches(:, 3) - Opt4(:, 1) ) .^ 2 );
mse_Opt(2,:) = mse_Opt(2,:) + sum( ( matches(:, 4) - Opt4(:, 2) ) .^ 2 );


