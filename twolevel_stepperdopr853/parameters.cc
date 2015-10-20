const Doub pi = acos(-1);
const int nvar = 3;	//	微分方程组未知数个数
static double omega_0;		//二能级原子跃迁频率,原子单位,5.36来自YangWeiFeng
static double omega_L = 0.056;	//基频场频率，原子单位 a.u.
static double rabbi_0;		//基频场拉比频率，原子单位
static double mu = 1.0e-29 / 1.60217653e-19 / 5.291772108e-11;		//偶极跃迁几率a.u.
static double xi;
static double rabbi_tmp;
static int laserchoice;
static double cycles_g;			//脉冲FWHM周期数
static double T;
static double dur;
static int IF_Chirp;
static double eta = 6.25;
static double tao = 120.0;
static double chirp_phase;
static int decaySwitch=1;
static double T1(1.0E-12 / 2.418884326505E-17);
static double T2(0.5E-12 / 2.418884326505E-17);
static double w_init;

static double Natoms=7.5E24*pow(5.29E-11, 3);
static int CollectFlag=1;