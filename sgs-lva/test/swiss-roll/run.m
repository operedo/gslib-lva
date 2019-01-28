
% based in scripts from http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.html

ppm = [4000 4000 4000 4000];        % points per mixture
const=1.0;
centers = const * [7.5 7.5; 7.5 12.5; 12.5 7.5; 12.5 12.5];
stdev = 3;
[data2,labels]=makegaussmixnd(centers,stdev,ppm);
%plotcol (data2, ppm, 'rgbk');

X = data2(:,1);
Y = 2*data2(:,2);

maxvalX=max(X);
minvalX=min(X);
maxvalY=max(Y);
minvalY=min(Y);

Xnorm = (X-minvalX*ones(size(X)))/(maxvalX-minvalX);
Ynorm = (Y-minvalY*ones(size(Y)))/(maxvalY-minvalY);

data2norm = [Xnorm, Ynorm];
plotcol (data2norm, ppm, 'rgbk');

turns=20
data3 = [Xnorm.*sin(Xnorm*pi*2*turns) Xnorm.*cos(Xnorm*pi*2*turns) Ynorm];
degrees = [rad2deg(Xnorm*pi*2*turns) rad2deg(Xnorm*pi*2*turns) zeros(size(X))];
plotcol (data3,ppm,'rgbk');

%data4 = [data3 degrees labels'];
data4 = [data3 degrees 10*ones(size(X)) ones(size(X))];
csvwrite('swiss_roll_v4.csv',data4);

data5 = [data3 labels'];
csvwrite('swiss_roll_v4_labels.csv',data5);
