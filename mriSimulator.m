%% Constants
larmor = 3*(42.6*(10^6));
h=6.6*10^(-34);
k=1.4*10^(-23);
per = 310.2;
t1map = zeros(217,181,181);
t2map = zeros(217,181,181);
pdmap = zeros(217,181,181);
t2map2 = zeros(217,181,181);
m0cons = ((pi*larmor)*(h)^2)/(4*(k*per));%m0 formula from the report

%Brainweb constants which we used for data elimination
%Constans for background
t1_0=0;
t2_0=0;
t2s_0=0;
pd_0=0;
%Constans for CSF
t1_1=2569;
t2_1=329;
t2s_1=58;
pd_1=1;
%Constans for grey matter
t1_2=833;
t2_2=83;
t2s_2=69;
pd_2=0.86;
%Constans for white matter
t1_3=500;
t2_3=70;
t2s_3=61;
pd_3=0.77;
%Constans for fat
t1_4=350;
t2_4=70;
t2s_4=58;
pd_4=1;
%Constans for muscle/skin
t1_5=900;
t2_5=47;
t2s_5=30;
pd_5=1.0;
%Constans for skin
t1_6=2569;
t2_6=329;
t2s_6=58;
pd_6=1;
%Constans for skull
t1_7=0;
t2_7=0;
t2s_7=0;
pd_7=0;
%Constans for glial matter
t1_8=833;
t2_8=83;
t2s_8=69;
pd_8=0.86;
%Constans for meat
t1_9=500;
t2_9=70;
t2s_9=61;
pd_9=0.77;
%Constans for MS lesion
t1_10=752;
t2_10=237;
t2s_10=204;
pd_10=0.76;

%Image dimensions
x=1:217;
y=1:181;
z=1:181;%Number of slices of brain

alpha=pi/2;%flip angle(TE=35 and TR=100)
b0 =zeros(217,181);%b0 inhomogenity

%For T1 weighted imaging: TE=30 and TR=1000
%For PD weighted imaing: TE=30 and TR=4000
%For T2 weighted imaging: TE=100 and TR=4000
TE=35;
TR=100;

f=zeros(217,181);
f3d=zeros(217,181,181);
f_fourier= zeros(217,181,181);
shiftedf= zeros(217,181,181);
newf= zeros(217,181,181);
%% Code
image = loadminc("phantom_1.0mm_normal_crisp.mnc");%%dataset from brainweb
% Data Elimination
for i = 1:217
    for m = 1:181
        for l = 1:181
            if (image(i,m,l)==0)%background
                t1map(i,m,l)=t1_0;
                t2map(i,m,l)=t2_0;
                t2map2(i,m,l)=t2s_0;
                pdmap(i,m,l)=pd_0;
            end
            if (image(i,m,l)==1)%CSF
                t1map(i,m,l)=t1_1;
                t2map(i,m,l)=t2_1;
                t2map2(i,m,l)=t2s_1;
                pdmap(i,m,l)=pd_1;
            end
            if (image(i,m,l)==2)%grey matter
                t1map(i,m,l)=t1_2;
                t2map(i,m,l)=t2_2;
                t2map2(i,m,l)=t2s_2;
                pdmap(i,m,l)=pd_2;
            end
            if (image(i,m,l)==3)%white matter              
                t1map(i,m,l)=t1_3;
                t2map(i,m,l)=t2_3;
                t2map2(i,m,l)=t2s_3;
                pdmap(i,m,l)=pd_3;
            end
            if (image(i,m,l)==4)%fat              
                t1map(i,m,l)=t1_4;
                t2map(i,m,l)=t2_4;
                t2map2(i,m,l)=t2s_4;
                pdmap(i,m,l)=pd_4;
            end
            if (image(i,m,l)==5)%muscle/skin                
                t1map(i,m,l)=t1_5;
                t2map(i,m,l)=t2_5;
                t2map2(i,m,l)=t2s_5;
                pdmap(i,m,l)=pd_5;
            end
            if (image(i,m,l)==6)%skin                
                t1map(i,m,l)=t1_6;
                t2map(i,m,l)=t2_6;
                t2map2(i,m,l)=t2s_6;
                pdmap(i,m,l)=pd_6;
            end
            if (image(i,m,l)==7)%skull               
                t1map(i,m,l)=t1_7;
                t2map(i,m,l)=t2_7;
                t2map2(i,m,l)=t2s_7;
                pdmap(i,m,l)=pd_7;
            end
            if (image(i,m,l)==8)%glial matter               
                t1map(i,m,l)=t1_8;
                t2map(i,m,l)=t2_8;
                t2map2(i,m,l)=t2s_8;
                pdmap(i,m,l)=pd_8;
            end
            if (image(i,m,l)==9)%meat                
                t1map(i,m,l)=t1_9;
                t2map(i,m,l)=t2_9;
                t2map2(i,m,l)=t2s_9;
                pdmap(i,m,l)=pd_9;
            end
            if (image(i,m,l)==10)%MS lesion               
                t1map(i,m,l)=t1_10;
                t2map(i,m,l)=t2_10;
                t2map2(i,m,l)=t2s_10;
                pdmap(i,m,l)=pd_10;
            end
        end
    end
end
M0=m0cons.*pdmap;
forg(x,y,z)=M0(x,y,z).*sin(alpha).*exp(-1.*TE./t2map(x,y,z)).*(1-exp(-1.*TR./t1map(x,y,z)))./(1-cos(alpha).*exp(-1.*TR./t1map(x,y,z)));
%Original Image
figure,imshow(forg(:,:,75),[]);
title("Original Image")
%Data Reconstruction
for i = 1:181
    f_fourier(:,:,i) = fft2im(forg(:,:,i));
    shiftedf(:,:,i) = ifftshift(ifftn(f_fourier(:,:,i)));
    newf(:,:,i) = squeeze(sqrt(sum(abs(shiftedf(:,:,i)).^2,42)));
end

figure,imshow(newf(:,:,75),[]);
title("Reconstructed Image")
%Tissue Contrast Values
%Make table of csf, grey_matter and white_matter values for TR=1000,10000 and TE=20,40,60,120,300,600,1200
%Make table of csf, grey_matter and white_matter values for TE=40,400 and TR=10000,5000,2500,1000,500,250,100
csf = newf(167,167,75);
grey_matter = newf(103,80,75);
white_matter = newf(77,152,75);

%B0 Inhomogenity(TE=35 and TR=100)
for i = 1:217
    for j = 1:181
        b0(i,j)=3*(1+0.01*sin(2*pi*i/217))*(1+0.01*cos(2*pi*j/181));
        b1(i,j)=3*(1+0.1*sin(2*pi*i/217))*(1+0.1*cos(2*pi*j/181));
    end
end
b_0=round(b0);
b_1=round(b1);
subplot(2,1,1);
f2(x,y)=b_0.*2.*pdmap(:,:,72).*sin(alpha).*exp(-1.*TE./t2map(x,y,72)).*(1-exp(-1.*TR./t1map(x,y,72)))./(1-cos(alpha).*exp(-1.*TR./t1map(x,y,72)));
imshow(f2);
title("B0 inhomogeneity1")
subplot(2,1,2);
f3(x,y)=b_1.*2.*pdmap(:,:,72).*sin(alpha).*exp(-1.*TE./t2map(x,y,72)).*(1-exp(-1.*TR./t1map(x,y,72)))./(1-cos(alpha).*exp(-1.*TR./t1map(x,y,72)));
imshow(f3);
title("B0 inhomogeneity2")
%% Functions

function res = fft2im(x)
    res = 1/sqrt(length(x(:)))*fftshift(fft2(ifftshift(x)));
end