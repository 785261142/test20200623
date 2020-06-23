%this is paper 4 code .
clc
clear all;
close all;

%%figure 2
% Monte Carlo simulation of two Exponential RVs
avg_y1=1;
avg_y2=2;
n_loop=100000;
MonteCarlo=[];
x=0.05:0.1:4-0.05;
y1=exprnd(avg_y1,[1 n_loop]);
y2=exprnd(avg_y2,[1 n_loop]);
y=y1.*y2./(y1+y2);
yy=y1.*y2./(y1+y2+1);
n=hist(y,x)./n_loop*10;
figure(1);
bar(x,n);
hold on
%贝塞尔函数
K1 = besselk(1,2*x./sqrt(avg_y1.*avg_y2));
K0 = besselk(0,2*x./sqrt(avg_y1.*avg_y2));
py=(2*x.*exp(-x.*(1/avg_y1+1/avg_y2))./(avg_y1.*avg_y2)).*((avg_y1+avg_y2)./sqrt(avg_y1.*avg_y2).*K1+2.*K0);
plot(x+0.05,py,'r');
legend('\itsimulation','\itanalysis');
ylabel('probability');
xlabel('SNR');

%%figure 3
% outage probability of two gains
yth=1;
for avg=1:40	%dB
	avg_y1=10^(avg/10);
	avg_y2=10^(avg/10);
	n_loop=1000000;
	y1=exprnd(avg_y1,[1 n_loop]);
	y2=exprnd(avg_y2,[1 n_loop]);
	y=y1.*y2./(y1+y2);
	yy=y1.*y2./(y1+y2+1);
	gain1(avg)=sum(y<yth)./n_loop;
	gain2(avg)=sum(yy<yth)./n_loop;
end 
figure(2),semilogy(gain1,'-*'),hold on;
semilogy(gain2,'-'),grid on;
legend('\itGain1','\itGain2');
ylabel('outage probability');
xlabel('SNR');

%%figure 4
yth=1;
for avg=1:40
	avg_y1=10^(avg/10);
	avg_y2=10^(avg/10);
	K1 = besselk(1,2*yth/sqrt(avg_y1*avg_y2));
	Poutnonregenerative(avg)=1-2*yth/sqrt(avg_y1*avg_y2)*K1*exp(-yth*(1/avg_y1+1/avg_y2));
	Poutregenerative(avg)=1-exp(-yth*(1/avg_y1+1/avg_y2));
end 
for avg=1:40
	avg_y1=20^(avg/10);
	avg_y2=10^(avg/10);
	K1 = besselk(1,2*yth/sqrt(avg_y1*avg_y2));
	Poutnonregenerative1(avg)=1-2*yth/sqrt(avg_y1*avg_y2)*K1*exp(-yth*(1/avg_y1+1/avg_y2));
	Poutregenerative2(avg)=1-exp(-yth*(1/avg_y1+1/avg_y2));
end 
figure(3),semilogy(Poutnonregenerative,'-r'),hold on;
semilogy(Poutregenerative,'-*g'),hold on
hold on
semilogy(Poutnonregenerative1,'-r'),hold on;
semilogy(Poutregenerative2,'-*g');
grid on;
legend('\itnonregenerative r1=r2','\itregenerative r1=r2','\itnonregenerative r1=2r2','\itregenerative r1=2r2');
ylabel('outage probability');
xlabel('SNR');
hold on
