
//两能级原子的高次谐波的数值程序，实际求解三元一阶常微分方程组（布洛赫方程）的初值问题，
//利用定步长积分四阶Runge-Kutta Method
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#include"../include/nr3.h"
#include"../include/odeint.h"
#include"../include/stepper.h"
#include"../include/stepperdopr853.h"
#include "parameters.cc"
#include "odeequations.cc"

namespace BFS = boost::filesystem;
using namespace std;
const int FloatPonitNumber=10;

template <class Ty>
inline void writeformat(Ty& os,int n=10)
{
	os << scientific << setprecision(n);
};

int main(int argc, char* argv[])
{
	const Doub atol = 1.0e-8, rtol = atol, hmin = 0.0;
	int i;
	//initial values//

	VecDoub ystart(nvar);
	cout << "分别输入布洛赫方程u、v和w的初始值：" << endl;
	for (i = 0; i < nvar; i++)
	{
		cin >> ystart[i];
	}
	w_init = ystart[nvar-1];

	///matter parameters///
	cout << "输入能级跃迁频率omega_0：";
	cin >> omega_0;

	double mu_11, mu_22;
	cout << "分别输入固有偶极矩mu_11/mu_22(单位以偶极跃迁矩阵元mu的倍数)： " << endl;
	cin >> mu_11;
	cin >> mu_22;
	mu_11 *= mu;
	mu_22 *= mu;
	xi = (mu_22 - mu_11) / (2.0*mu);

	//////Laser parameters /////
	cout << "输入驱动脉冲拉比频率rabbi= ";
	cin >> rabbi_0;
	T = 2 * pi / omega_L;			//光周期， 单位 a.u.
	cout << "矩形场输入数字0，高斯场输入数字1，sin-square输入2：";
	cin >> laserchoice;
	cout << "输入脉冲的FWHM周期数cycles_g：";
	cin >> cycles_g;
	dur = cycles_g * T;
	double cycles;				//脉冲光学周期数
	if (laserchoice == 0)	//矩形
	{
		cycles = cycles_g + 2.0;
	}
	else if (laserchoice == 1)	//高斯
	{
		cycles = cycles_g * 3;
	}
	else if (laserchoice == 2)	//sin-square
	{
		cycles = cycles_g*2.0 + 2.0;
	}

	int n = 13;		//
	int N = pow(2, n);				//每个光周期的积分步长数
	int NN = N*cycles;				//整个脉冲积分总步长数，不包括起始点
	int NT = NN + 1;					//整个脉冲积分总步长数，包括起始点
	double h = T / double(N);			//积分步长
	double half_h = h / 2.0;
	double tstart = -cycles / 2.0 * T;
	double tend = -tstart;
	bool InputAgainFlag=false;
	do 
	{
		cout << "是否有啁啾，数字1有，数字0没有：";
		cin >> IF_Chirp;
		if (IF_Chirp != 0 && IF_Chirp != 1)
		{
			InputAgainFlag=true;
		}
		else
		{
			InputAgainFlag=false;
		}
	} while (InputAgainFlag==true);
	
	do 
	{
		cout << "偶极矩的计算是否考虑集体效应，数字1考虑，数字0不考虑： ";
		cin >> CollectFlag;
		if (CollectFlag != 0 && CollectFlag != 1)
		{
			InputAgainFlag=true;
		}
		else
		{
			InputAgainFlag=false;
		}
	} while (InputAgainFlag==true);	
	if (CollectFlag == 0)
	{
		Natoms=1;
	}

	do 
	{
		cout << "是否考虑衰减，数字1考虑，数字0不考虑： ";
		cin >> decaySwitch;
		if (decaySwitch != 0 && decaySwitch != 1)
		{
			InputAgainFlag=true;
		}
		else
		{
			InputAgainFlag=false;
		}
	} while (InputAgainFlag==true);
	

	///// Laser field /////	
	VecDoub laser_field(NT);								//记录脉冲对应t时刻的电场值
	VecDoub rabbi(NT);
	VecDoub rabbi_deriv(NT);

	BFS::path current_dir = BFS::current_path();
	BFS::path path_res = current_dir.parent_path() / "res";
	if (!BFS::exists(path_res))
	{
		cout << "There's no directory " << path_res.filename() << ". Program will create it." << endl;
		BFS::create_directory(path_res);
	}

	BFS::ofstream time_output(path_res / "time.txt", std::ofstream::out | std::ofstream::trunc);
	BFS::ofstream lasersource_output(path_res / "source.txt", std::ofstream::out | std::ofstream::trunc);

	writeformat(time_output,FloatPonitNumber);
	writeformat(lasersource_output,FloatPonitNumber);

	double t = tstart;						//表征脉冲时刻
	for (i = 0; i < (NT - 1); i++)
	{
		if (IF_Chirp == 0)
		{
			chirp_phase = 0;
		}
		else
		{
			chirp_phase = -eta*tanh(t / tao);
		}

		if (laserchoice == 0)
		{
			if (t >= -cycles_g / 2.0*T&&t <= cycles_g / 2.0*T)
			{
				laser_field[i] = rabbi_0 / mu*sin(omega_L*t);
			}
			else
			{
				laser_field[i] = 0.0;
			}
		}
		else if (laserchoice == 1)
		{
			laser_field[i] = rabbi_0 / mu*exp(-4 * log(2)*t*t / dur / dur)*cos(omega_L*t + chirp_phase);	//高斯型激光场
		}
		else if (laserchoice == 2)
		{
			if (t >= -cycles_g*T&&t <= cycles_g*T)
			{
				laser_field[i] = rabbi_0 / mu*pow(sin(pi*(t + cycles_g*T) / (cycles_g*2.0*T)), 2)*cos(omega_0*t + chirp_phase);
			}
			else
			{
				laser_field[i] = 0.0;
			}
		}

		rabbi[i] = mu*laser_field[i];
		lasersource_output << laser_field[i] << endl;

		t = t + h;
	}
	lasersource_output.close();

	//equation propagate//
	Output out(NN); //Dense output at NT-1 points plus x1.
	rhs_van d; //Declare d as a rhs_van object.
	Odeint<StepperDopr853<rhs_van> > ode(ystart, tstart, tend, atol, rtol, h, hmin, out, d);
	ode.integrate();

	//将相关参数写入文件
	BFS::fstream relative_parameters;

	relative_parameters.open(path_res / "relative_parameters.dat", std::ios_base::out | std::ios_base::binary);
	relative_parameters.write((char*)&omega_0, sizeof(double));
	relative_parameters.write((char*)&omega_L, sizeof(double));
	relative_parameters.write((char*)&rabbi_0, sizeof(double));
	relative_parameters.write((char*)&mu, sizeof(double));
	relative_parameters.write((char*)&xi, sizeof(double));
	relative_parameters.write((char*)&NT, sizeof(int));
	relative_parameters.write((char*)&h, sizeof(double));
	relative_parameters.write((char*)&T, sizeof(double));
	relative_parameters.write((char*)&tstart, sizeof(double));
	relative_parameters.close();


	/////将平均偶极矩计算结果写入文件dipole.dat/////	
	BFS::ofstream dipole(path_res / "dipole.txt", std::ofstream::out | std::ofstream::trunc);
	writeformat(dipole,FloatPonitNumber);
	for (i = 0; i < out.count - 2; i++)
	{
		time_output << out.xsave[i] / T<<endl;
		dipole << Natoms*(mu*out.ysave[0][i] + mu_11*(1.0 - out.ysave[2][i])/ 2.0 + 
					mu_22*(1.0 + out.ysave[2][i]) / 2.0)
					<< " " << Natoms*mu*out.ysave[0][i] <<endl;
	}
	time_output.close();
	dipole.close();

	return 0;
}

