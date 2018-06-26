#ifndef aes_h_
#define aes_h_

typedef unsigned char byte;

static const byte sbox[256];
static const byte rsbox[256];
static const byte Rcon[10];
static const byte M2[256];
static const byte M3[256];
static const byte M9[256];
static const byte M11[256];
static const byte M13[256];
static const byte M14[256];

byte *XOR_matrix(byte M1[16],byte M2[16]);
byte Mult2(byte a);
byte Mult3(byte a);
byte Mult9(byte a);
byte Mult11(byte a);
byte Mult13(byte a);
byte Mult14(byte a);
byte get_sbox(byte b);
byte get_rsbox(byte b);

void printState(byte X[16]);
void print2States(byte X[16],byte Y[16]);
int equals(byte C[16],byte C2[16]);

void KeySchedule(byte K[16],int round);
void KeyExpansion(byte K[16],byte RK[11][16]);
void KeySchedule_1(byte K[16],int round);
void AddRoundKey(byte X[16], byte RK[11][16],int round);
void AddRoundKey_1(byte X[16], byte K[16],int round);
void SubBytes(byte X[16]);
void SubBytes_1(byte X[16]);
void ShiftRows(byte X[16]);
void ShiftRows_1(byte X[16]);
void MixColumns(byte X[16]);
void MixColumns_1(byte X[16]);

void AES_128(byte X[16],byte K[16],byte Y[16]);
void AES_128_1(byte X[16],byte K[16],byte Y[16]);

void BicliqueKeyExpansion(byte K[16],byte RK[3][16]);
void BicliqueKeyExpansion_1(byte K[16],byte RK[9][16]);
void g(byte X[16],byte K[16],byte Y[16]);
void f(byte X[16],byte K[16],byte Y[16]);
void f_1(byte X[16],byte K[16],byte Y[16]);
byte r(byte X[16],byte precomp[7],byte RK[3][16]);
byte t_1(byte X[16],byte precomp[2][16],byte RK[9][16]);

#endif  /*aes_h*/
