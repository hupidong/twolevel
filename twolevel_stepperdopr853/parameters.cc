const Doub pi = acos(-1);
const int nvar = 3;	//	΢�ַ�����δ֪������
static double omega_0;		//���ܼ�ԭ��ԾǨƵ��,ԭ�ӵ�λ,5.36����YangWeiFeng
static double omega_L = 0.056;	//��Ƶ��Ƶ�ʣ�ԭ�ӵ�λ a.u.
static double rabbi_0;		//��Ƶ������Ƶ�ʣ�ԭ�ӵ�λ
static double mu = 1.0e-29 / 1.60217653e-19 / 5.291772108e-11;		//ż��ԾǨ����a.u.
static double xi;
static double rabbi_tmp;
static int laserchoice;
static double cycles_g;			//����FWHM������
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