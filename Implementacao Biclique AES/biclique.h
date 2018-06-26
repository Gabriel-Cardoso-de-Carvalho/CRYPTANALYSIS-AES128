#ifndef biclique_h_
#define biclique_h_

typedef unsigned char byte;

void forward_keys(byte RK[9][16],byte RK_ret[3][16], int j);
byte forward_r(byte X[16],byte fullforward[7],byte RK[3][16]);
void backward_keys(byte RK[9][16],byte RK_ret[5][16], int i);
byte backward_t_1(byte X[16],byte fullbackward[2][16],byte RK[5][16]);

void ConstructBiclique(byte C0[16],byte K0[16],byte Ci[256][16],byte Sj[256][16]);
void GeneratePi(byte real_K[16],byte Pi[256][16],byte Ci[256][16]);
int MITM(byte K[16],byte Pi[256][16],byte Sj[256][16],byte real_K[16]);
int MP(byte K[16],byte Pi[256][16],byte Sj[256][16],byte real_K[16]);
void Precomputation(byte K[16],byte Pi[256][16],byte Sj[256][16]);
int Recomputation(byte K[16],byte Pi[256][16],byte Sj[256][16],byte res[16]);
int BicliqueAttack_MITM(byte C0[16],byte real_K[16], byte K0[16]);
int BicliqueAttack_MP(byte C0[16],byte real_K[16], byte K0[16]);


#endif  /*biclique_h*/
