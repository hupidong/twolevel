
struct rhs_van
{
	//	rhs_van(Int lc) : laserchoice(lc) {}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
	{
		chirp_phase = -IF_Chirp*eta*tanh(x / tao);
		if (laserchoice == 0)
		{
			if (x >= -cycles_g / 2.0*T&&x <= cycles_g / 2.0*T)
			{
				rabbi_tmp = rabbi_0*sin(omega_L*x);
			}
			else
			{
				rabbi_tmp = 0.0;
			}
		}
		else if (laserchoice == 1)
		{
			rabbi_tmp = rabbi_0*exp(-4.0 * log(2.0)*x*x / dur / dur)*cos(omega_L*x + chirp_phase);	//高斯型激光场
		}
		else if (laserchoice == 2)
		{
			if (x >= -cycles_g*T&&x <= cycles_g*T)
			{
				rabbi_tmp = rabbi_0*pow(sin(pi*(x + cycles_g*T) / (cycles_g*2.0*T)), 2)*cos(omega_L*x + chirp_phase);
			}
			else
			{
				rabbi_tmp = 0.0;
			}
		}
		if (decaySwitch == 0)
		{
			dydx[0] = -omega_0*y[1] + 2.0*xi*rabbi_tmp * y[1];
			dydx[1] = omega_0*y[0] - 2.0*xi*rabbi_tmp * y[0] + 2.0*rabbi_tmp * y[2];
			dydx[2] = -2.0*rabbi_tmp * y[1];
		}
		else if (decaySwitch == 1)
		{
			dydx[0] = -omega_0*y[1] + 2.0*xi*rabbi_tmp * y[1]-y[0]/T2;
			dydx[1] = omega_0*y[0] - 2.0*xi*rabbi_tmp * y[0] + 2.0*rabbi_tmp * y[2]-y[1]/T2;
			dydx[2] = -2.0*rabbi_tmp * y[1]-(y[2]-w_init)/T1;
		}
		else
		{
			throw("是否考虑衰减的选择输入有错误！！！");
		}
	}
};