#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <complex>

namespace BFS=boost::filesystem;
const double pi = acos(-1.0);

int main(int argc, char*argv[])
{
	using namespace std;
	using VecDoub = std::vector<double>;
	using MatDoub = std::vector < std::vector<double> >;
	using VecComplex = std::vector < std::complex<double> >;
	using MatComplex = std::vector < std::vector<std::complex<double>>>;
	int i = 0, j = 0, k = 0;

	BFS::path current_dir=BFS::current_path();
	BFS::path path_res = current_dir.parent_path() / "res";
	BFS::path current_file;
	if (!BFS::exists(path_res))
	{
		cout << "There's no directory " << path_res.filename() << ". Program will be terminated." << endl;
		system("pause");
		exit(0);
	}

	double omega_0, omega_L, rabbi_0, mu, T, tao, tstart, dt,xi;
	int N;

	current_file = path_res / "relative_parameters.dat";
	BFS::fstream para;
	para.open(current_file, ios::in | ios::binary);
	if (!para)
	{
		cout << "Can not open file " << current_file.string() << endl;
		system("pause");
		exit(0);
	}
	para.read((char *)&omega_0, sizeof(double));
	para.read((char *)&omega_L, sizeof(double));
	para.read((char *)&rabbi_0, sizeof(double));
	para.read((char *)&mu, sizeof(double));
	para.read((char *)&xi, sizeof(double));
	para.read((char*)&N, sizeof(int));
	para.read((char*)&dt, sizeof(double));
	para.read((char*)&T, sizeof(double));
	para.read((char *)&tstart, sizeof(double));
	para.close();
	N -= 1;
	tao = 6.0;//wavelet parameter

	VecDoub t(N);
	for (i = 0; i<N; i++)
	{
		t[i] = tstart + i*dt;
	}

//	double wti, wtf, wdt, wwi, wwf, wdw;//小波分析的起始时间，结束时间，步长，最小频率，最大频率，频率间隔
	int reSampleScale = pow(2, 3);
	int sampleNum(0);
	VecDoub tSample;
	for (i = 0; i < N; i+=reSampleScale)
	{
		tSample.push_back(t[i]);
		sampleNum++;
	}
	
	double tSelectBegin(tstart), tSelectEnd, Dt(reSampleScale*dt);
	double wLower,wUpper,Dw;
	double scaleLower,scaleUpper,dScale;
	/*
	cout << "The total time start from " << tstart / T << " T, and end at time " << t[N - 1] / T << " T."<<endl;
	cout << "Select the begin time point for wavelet (unit in T): ";
	cin>>tSelectBegin;
	tSelectBegin*=T;
	cout << "Then select the end time point for wavelet (unit in T): ";
	cin>>tSelectEnd;
	tSelectEnd*=T;
	*/
	cout << "Now input the lower harmonic order for wavelet (unit in order): ";
	cin>>wLower;
	wLower*=omega_L;
	cout << "Then input the upper harmonic order for wavelet (unit in order): ";
	cin>>wUpper;
	wUpper*=omega_L;
	int scaleNum;
	cout << "And, input the scale numbers for wavelet: ";
	cin>>scaleNum;
	scaleUpper=1.0/wLower;
	scaleLower=1.0/wUpper;

	int nt(sampleNum), nw(scaleNum);
	cout << "nt=" << nt << endl; cout << "nw=" << nw << endl;
	VecDoub tcenter(nt), freqs(nw), scales(nw);
	current_file = path_res / "wave_tcenter.txt";
	BFS::ofstream wave_tcenter(current_file, std::ofstream::out | std::ofstream::trunc);
	current_file = path_res / "wave_freqs.txt";
	BFS::ofstream wave_freqs(current_file, std::ofstream::out | std::ofstream::trunc);
	double unit = 0;

	tcenter=tSample;
	for (i = 0; i<nt; i++)
	{
		unit = tcenter[i] / T;
		wave_tcenter<<unit<<endl;
	}
	dScale = (scaleUpper-scaleLower) / (nw-1);
	for (i = 0; i<nw; i++)
	{
		scales[i] = scaleUpper-i*dScale;
		freqs[i] = 1.0 / scales[i];
		unit = freqs[i] / omega_L;
		wave_freqs<<unit<<endl;
	}
	wave_tcenter.close();
	wave_freqs.close();

	MatComplex Coeffs(nw, VecComplex(nt));
	VecDoub dipole(N), dipoleSample(nt);

	BFS::ifstream dipoleInput;
	current_file = path_res / "dipole.txt";

	dipoleInput.open(current_file, std::ifstream::in);
	if (!dipoleInput)
	{
		cout << "Can not open file " << current_file.string()<< endl;
		system("pause");
		exit(0);
	}
	double tmp1,tmp2;
	for (i = 0; i<N; i++)
	{
		dipoleInput >> tmp1 >> tmp2;
		dipole[i]=tmp1;
	}
	for (i = 0,j=0; i < N; i+=reSampleScale,j++)
	{
		dipoleSample[j] = dipole[i];
	}
	dipoleInput.close();

	/*	Engine *ep;
	mxArray *tt=NULL,*result=NULL,*ts;
	if(!(ep=engOpen("\0")))
	{
	cout<<"Can not start matlab engine"<<endl;
	}
	tt=mxCreateDoubleMatrix(1,N,mxREAL);
	ts=mxCreateDoubleMatrix(1,N,mxREAL);
	memcpy((void *) mxGetPr(tt),(void*)datime,N*sizeof(double));
	engPutVariable(ep,"tt",tt);
	memcpy((void *) mxGetPr(ts),(void*)da,N*sizeof(double));
	engPutVariable(ep,"ts",ts);
	engEvalString(ep,"plot(tt,ts);");
	*/
	BFS::ofstream Coeffs_Realout,Coeffs_Imagout;
	cout << "ok" << endl;
	current_file = path_res / "Coeffs_Real.txt";
	Coeffs_Realout.open(current_file,std::ofstream::out|std::ofstream::trunc);
	if (!Coeffs_Realout)
	{
		cout << "Can not open file " << current_file.string() << endl;
		system("pause");
		exit(0);
	}
	current_file = path_res / "Coeffs_Imag.txt";
	Coeffs_Imagout.open(current_file, std::ofstream::out | std::ofstream::trunc);
	if (!Coeffs_Imagout)
	{
		cout << "Can not open file " << current_file.string() << endl;
		system("pause");
		exit(0);
	}
	double omgt;
	complex<double> I_Complex(0,1);
	for (i = 0; i<nw; i++)
	{
		for (j = 0; j<nt; j++)
		{
			Coeffs[i][j] = 0.0;
			for (k = 0; k<nt; k++)
			{
				omgt = freqs[i] * (tSample[k] - tcenter[j]);
				Coeffs[i][j] += dipoleSample[k] * sqrt(freqs[i] / tao)*exp(I_Complex*omgt)*exp(-omgt*omgt / 2.0 / tao / tao)*Dt;
			}
		}
	}
	for (i = 0; i < nw; i++)
	{
		for (j = 0; j < nt; j++)
		{
			if (j == (nt - 1))
			{
				Coeffs_Realout << Coeffs[i][j].real();
				Coeffs_Imagout << Coeffs[i][j].imag();
			}
			else
			{
				Coeffs_Realout << Coeffs[i][j].real() << " ";
				Coeffs_Imagout << Coeffs[i][j].imag() << " ";
			}
		}
		Coeffs_Realout << endl;
		Coeffs_Imagout << endl;
	}
	Coeffs_Realout.close();
	Coeffs_Imagout.close();

return 0;
}