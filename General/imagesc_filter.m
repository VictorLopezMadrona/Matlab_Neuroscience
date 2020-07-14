function IM_f=imagesc_filter(im_in,smooth,x_axis,y_axis)

%
% Plot a smooth image increasing the number of pixels and filtering with a 
% low-pass.
%
% USE:
%    imagesc_filter(im_in,smooth,x,y)
%
% INPUT:
%    im_in: matrix (MxN) with the image to plot.
%    smooth: Integer with the intensity of the smooth (Default = 0 -no smooth).
%    x,y (Optional): picture axis with format x=[x_min x_max].

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jun. 2016; Last revision: 14-Jul-2020

if nargin ==1
    n = 0;
elseif nargin==2
    n = smooth;
    [X,Y]=size(im_in);
    y_axis=[0 X];
    x_axis=[0 Y];
elseif nargin>2
    n = smooth;
end

if n==0
    imagesc(im_in)
    return
end

h=1/9*[1 1 1; 1 1 1; 1 1 1];
M_pre=max(max(im_in));
im_in=imfilter(imfilter(im_in,h),h);
M_post=max(max(im_in));
im_in=im_in*M_pre/M_post;

[X,Y]=size(im_in);
im_out=zeros(X*n,Y*n);

for x=1:X
    for y=1:Y
        im_out((x-1)*n+1:x*n, (y-1)*n+1:y*n) = im_in(x,y);
    end
end

IM_f=imfilter(imfilter(im_out,h),h);

xx=linspace(x_axis(1),x_axis(2),X*4);
yy=linspace(y_axis(1),y_axis(2),Y*4);
imagesc(xx,yy,IM_f); axis('xy')

        