pkg load signal
u = double(imread('../images/lena.png'));
H = hist(u(:),0:255);
idx=FTC_Seg(H,0);
