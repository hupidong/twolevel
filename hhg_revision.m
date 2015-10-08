clear;
data_path='C:\Users\PeterHu-T420\Desktop\HP_DESKTOP\twolevel\res';
parameter_path=strcat(data_path,'\relative_parameters.dat');

fid=fopen(parameter_path,'r');
omega0=fread(fid,1,'double');
omegaL=fread(fid,1,'double');
Omega0=fread(fid,1,'double');
mu=fread(fid,1,'double');
xi=fread(fid,1,'double');
len=fread(fid,1,'int32');
dt=fread(fid,1,'double'); %a.u.
T=fread(fid,1,'double');
fclose(fid);
wmg=2*pi/dt;

%plot laser field && fft 
formatDouble='%f64';
source_path=strcat(data_path,'\source.txt');
source_fread=fopen(source_path,'r');
source=textscan(source_fread,formatDouble);
source=cell2mat(source);

dipole_path=strcat(data_path,'\dipole.txt');
dipole_fread=fopen(dipole_path,'r');
dipole=textscan(dipole_fread,'%f64 %f64');
dipole=cell2mat(dipole);
dipole=dipole(:,1);
fclose(source_fread);fclose(dipole_fread);

time_path=strcat(data_path,'\time.txt');
time_fread=fopen(time_path,'r');
t=textscan(time_fread,formatDouble);
t=cell2mat(t);
fclose(time_fread);

figure
plot(t,source,'b-','linewidth',2);
title('Laser field');
xlabel('Time(In Optical Period )','fontsize',14);
ylabel('\it{E(t)}(a.u.)','fontsize',14);

%fft && plot
% if collection_effect==1
%     dipole=dipole*7.5E24*(5.29E-11)^3;
% end
len2=length(dipole);
figure
plot(t,dipole,'r-','linewidth',2)
title('dipole');
xlabel('Time(In Optical Period)','fontsize',14);
ylabel('Dipole Intensity (a.u.)','fontsize',14);

fre=(0:round(len2/2)-1)/len2*wmg/omegaL;
FFA=fft(dipole);
hff=FFA;
FFA=abs(FFA(1:round(len2/2))/length(dipole));%
FFA=2*log10(abs(FFA)+eps);
figure
plot(fre,FFA,'r-','linewidth',2);
title('HHG');
xlabel('Harmonic Order(\omega/\omega_L)','fontsize',14);
ylabel('Harmonic Intensity(arb.unit)','fontsize',14);

%wavelet transform using matlab
%考虑尺度如何采样的问题，尺度平均的话，频率就不平均，反之亦然
%%resample first and freq units transform
index=1:2^2:len2;
tSample=t(index);   %单位，T
% tSample=tSample*T*2.418884326505E-17;   %单位，s
dipoleSample=dipole(index);%a.u. 
Fs=2^13/2^2*(3/8*1.0E15);    %采样频率,考虑重采样
fc=centfrq('cmor1-1');   %Hz
freqL=input('Input the Lower limit of freqrange (unit in order): ');
freqU=input('Input the Upper limit of freqrange (unit in order): ');
freqL=freqL*omegaL*(1.0/2.418884326505E-17)/(2*pi);  %a.u. to Hz
freqU=freqU*omegaL*(1.0/2.418884326505E-17)/(2*pi);  %a.u. to Hz
freqrange=[freqL freqU];
scalerange=fc./(freqrange*(1/Fs));
scalesNum=256; %尺度采样个数
freqs=linspace(freqL,freqU,scalesNum);
scales=fc./(freqs*(1/Fs));
%scales=linspace(scalerange(end),scalerange(1),scalesNum);
Coeffs=cwt(dipoleSample,scales,'cmor1-1');
figure;
SCImg = wscalogram('image',abs(Coeffs),'scales',scales,'ydata',dipoleSample,'xdata',tSample);

%freqs=scal2frq(scales,'cmor1-1',1/Fs)*(2*pi)*2.418884326505E-17/omegaL;
freqs=freqs*(2*pi)*2.418884326505E-17/omegaL;
[Freqs,Tcenter]=meshgrid(freqs,tSample);
%%%%%
harmonic_energy=sqrt((2.0*xi*mu*source-omega0).^2+4.0*(mu*source).^2)/omegaL;
harmonic_energy=harmonic_energy(index);
figure;
subplot(3,1,1);
plot(tSample,harmonic_energy);
subplot(3,1,2);
surf(Tcenter,Freqs,abs(Coeffs)');shading('interp');view(0,90);

%%%%%%%%%%%%%%%%%%

%%wavelet transform using C++
freqs_path=strcat(data_path,'\wave_freqs.txt');
tcenter_path=strcat(data_path,'\wave_tcenter.txt');
freqs_fread=fopen(freqs_path,'r');
tcenter_fread=fopen(tcenter_path,'r');
freqs=textscan(freqs_fread,formatDouble);
freqs=cell2mat(freqs);
tcenter=textscan(tcenter_fread,formatDouble);
tcenter=cell2mat(tcenter);
nw=length(freqs);nt=length(tcenter);
fclose(freqs_fread);fclose(tcenter_fread);

Coeffs_Real_path=strcat(data_path,'\Coeffs_Real.txt');
Coeffs_Real_fread=fopen(Coeffs_Real_path,'r');
Coeffs_Real=textscan(Coeffs_Real_fread,repmat('%f64',[1,nt]));
Coeffs_Real=cell2mat(Coeffs_Real);
Coeffs_Imag_path=strcat(data_path,'\Coeffs_Imag.txt');
Coeffs_Imag_fread=fopen(Coeffs_Imag_path,'r');
Coeffs_Imag=textscan(Coeffs_Imag_fread,repmat('%f64',[1,nt]));
Coeffs_Imag=cell2mat(Coeffs_Imag);
fclose(Coeffs_Real_fread);fclose(Coeffs_Imag_fread);
Coeffs=Coeffs_Real+Coeffs_Imag*1i;
[Freqs,Tcenter]=meshgrid(freqs,tcenter);
% figure;
subplot(3,1,3)
surf(Tcenter,Freqs,abs(Coeffs)');shading('interp');view(0,90);