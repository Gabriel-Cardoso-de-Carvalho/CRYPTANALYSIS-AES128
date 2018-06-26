#include "aes.h"
#include "biclique.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "string.h"
#include "math.h"

byte dummy[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
byte delta_k_i[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};///bytes 8 e 12
byte nabla_k_j[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};///bytes 1 e 9
byte vi[256][256],vj,fullforward[256][7],fullbackward[256][4][16], RKi[256][9][16], RKj[256][9][16];

int is_candidate(byte v,int i,int j){
	return v == vi[i][j];
}

int check_candidate(byte P[16], byte K[16], byte S[16], int i, int j){
	byte testS[16],testK[16];
	
	memcpy(testK,K,16);
	
	testK[1] ^= j;
	testK[9] ^= j;
	testK[8] ^= i;
	testK[12]^= i;
	
	g(P,testK,testS);
	
	return equals(S,testS);
}

void forward_keys(byte RK[9][16],byte RK_ret[3][16], int j){
	byte a,b,c,d,A,B,C,D,E;
	
	memcpy(RK_ret[0],RK[0],16);
	memcpy(RK_ret[1],RK[1],16);
	memcpy(RK_ret[2],RK[2],16);
	
	a = (get_sbox(RK[7][13] ^ j) ^ RK[8][0] ^ 0x80) ^ RK[7][0];	
	b = (get_sbox(RK[4][12] ^ a) ^ RK[5][3]) ^ RK[4][3];
	c = (get_sbox(RK[3][13] ^ j) ^ (RK[4][0] ^ a) ^ 0x08) ^ RK[3][0];	
	D = get_sbox(RK[1][15] ^ b) ^ RK[2][2];
	d = D ^ RK[1][2];
	E = get_sbox(RK[0][12] ^ c) ^ (RK[1][3] ^ b);

	RK_ret[1][0] ^= c;
	RK_ret[1][1] ^= j;
	RK_ret[1][2]  = D;
	RK_ret[1][3] ^= b;
	RK_ret[1][5] ^= j;
	RK_ret[1][7] ^= b;
	RK_ret[1][8] ^= c;
	RK_ret[1][11]^= b;
	RK_ret[1][15]^= b;
	
	RK_ret[0][0] ^= c;
	RK_ret[0][1] ^= j;
	RK_ret[0][2] ^= d;
	RK_ret[0][3]  = E;
	RK_ret[0][4] ^= c;
	RK_ret[0][6] ^= d;
	RK_ret[0][8] ^= c;
	RK_ret[0][9] ^= j;
	RK_ret[0][12]^= c;
	
}

byte forward_r(byte X[16],byte fullforward[7],byte RK[3][16]){
	int i,j=0;
	byte Y[16];

	memcpy(Y,X,16);

	///AddRoundKey(Y,RK,0);
	Y[0] ^= RK[0][0];
	Y[1] ^= RK[0][1];
	Y[2] ^= RK[0][2];
	Y[3] ^= RK[0][3];
	Y[4] ^= RK[0][4];
	Y[6] ^= RK[0][6];
	Y[8] ^= RK[0][8];
	Y[9] ^= RK[0][9];
	Y[12]^= RK[0][12];

	/// SubBytes(Y);	
	Y[0] = get_sbox(Y[0]);
	Y[1] = get_sbox(Y[1]);
	Y[2] = get_sbox(Y[2]);
	Y[3] = get_sbox(Y[3]);
	Y[4] = get_sbox(Y[4]);
	Y[6] = get_sbox(Y[6]);
	Y[8] = get_sbox(Y[8]);
	Y[9] = get_sbox(Y[9]);
	Y[12]= get_sbox(Y[12]);
	
	/// ShiftRows(Y); Instead, read de value	
	Y[5] = Y[9];
	Y[7] = Y[3];
	Y[10] = Y[2];
	Y[13] = Y[1];
	Y[14] = Y[6];
	Y[1] = fullforward[0];
	Y[2] = fullforward[1];
	Y[3] = fullforward[2];
	Y[6] = fullforward[3];
	Y[9] = fullforward[4];
	Y[11]= fullforward[5];
	Y[15]= fullforward[6];
	
	MixColumns(Y);	
	
	/// AddRoundKey(Y,RK,1);
	Y[1] ^= RK[1][1];
	Y[6] ^= RK[1][6];
	Y[11]^= RK[1][11];
	Y[12]^= RK[1][12];
	

	/// SubBytes(Y);	
	Y[1] = get_sbox(Y[1]);
	Y[6] = get_sbox(Y[6]);
	Y[11]= get_sbox(Y[11]);
	Y[12]= get_sbox(Y[12]);
	
	/// ShiftRows(Y);	
	Y[13] = Y[1];
	Y[14] = Y[6];
	Y[15] = Y[11];
	
	/// MixColumns(Y);AddRoundKey(Y,RK,0);				
	return (Mult2(Y[12])^Mult3(Y[13])^Y[14]^Y[15]) ^ RK[2][12];
	
}

void backward_keys(byte RK[9][16],byte RK_ret[5][16], int i){

	memcpy(RK_ret[0],RK[3],16);
	memcpy(RK_ret[1],RK[4],16);
	memcpy(RK_ret[2],RK[5],16);
	memcpy(RK_ret[3],RK[6],16);
	memcpy(RK_ret[4],RK[7],16);
	
	byte a,A,B;
	

	
	A = get_sbox(RK[6][12] ^i) ^ RK[7][3];
	a = A ^ RK[6][3];
	
	B = get_sbox(RK[4][12] ^ i) ^ RK[5][3] ^ a;
	
	RK_ret[4][8] ^= i;

	RK_ret[3][3] = A;
	RK_ret[3][8] ^= i;
	RK_ret[3][12] ^= i;

	RK_ret[2][3] ^= a;
	RK_ret[2][8] ^= i;
	RK_ret[2][7] ^= a;

	RK_ret[1][3] = B;	
	RK_ret[1][8] ^= i;
	RK_ret[1][11]^= a;
	RK_ret[1][12]^= i;

	RK_ret[0][15]^= a;	
}

byte backward_t_1(byte X[16],byte fullbackward[4][16],byte RK[5][16]){
	int i,j=0,d[4];
	byte Y[16];
	
	memcpy(Y,X,16);

	AddRoundKey(Y,RK,4);
	
	///MixColumns_1(Y);
	memcpy(Y,fullbackward[0],8);
	memcpy(Y+12,fullbackward[0]+12,4);
	d[0] = Mult14(Y[8])^Mult11(Y[9])^Mult13(Y[10])^ Mult9(Y[11]);
	d[1] =  Mult9(Y[8])^Mult14(Y[9])^Mult11(Y[10])^Mult13(Y[11]);
	d[2] = Mult13(Y[8])^ Mult9(Y[9])^Mult14(Y[10])^Mult11(Y[11]);
	d[3] = Mult11(Y[8])^Mult13(Y[9])^ Mult9(Y[10])^Mult14(Y[11]);
	
	Y[8] = d[0];
	Y[9] = d[1];
	Y[10]= d[2];
	Y[11]= d[3];
		
	/// ShiftRows_1(Y);	
	Y[2] = Y[10];
	Y[7] = Y[11];
	Y[13]= Y[9];	
	memcpy(Y,fullbackward[1],2);
	memcpy(Y+3,fullbackward[1]+3,4);
	memcpy(Y+9,fullbackward[1]+9,4);
	memcpy(Y+14,fullbackward[1]+14,2);
	
	/// SubBytes_1(Y);	
	Y[2] = get_rsbox(Y[2]);
	Y[7] = get_rsbox(Y[7]);
	Y[8] = get_rsbox(Y[8]);
	Y[13]= get_rsbox(Y[13]);
	memcpy(Y,fullbackward[2],2);
	memcpy(Y+3,fullbackward[2]+3,4);
	memcpy(Y+9,fullbackward[2]+9,4);
	memcpy(Y+14,fullbackward[2]+14,2);
		
	/// AddRoundKey(Y,RK,3);
	Y[2] ^= RK[3][2];
	Y[3] ^= RK[3][3];
	Y[7] ^= RK[3][7];
	Y[8] ^= RK[3][8];
	Y[12]^= RK[3][12];
	Y[13]^= RK[3][13];
	memcpy(Y,fullbackward[3],2);
	memcpy(Y+4,fullbackward[3]+4,3);
	memcpy(Y+9,fullbackward[3]+9,3);
	memcpy(Y+14,fullbackward[3]+14,2);
	
	
	for(i = 2; i >=1;i--){
		
		MixColumns_1(Y);
		
		ShiftRows_1(Y);	
		
		SubBytes_1(Y);	
		
		AddRoundKey(Y,RK,i);
	}
	
	///Penultima
	MixColumns_1(Y);
		
	ShiftRows_1(Y);	
		
	/// SubBytes_1(Y);	
	
	Y[12] = get_rsbox(Y[12]);
	Y[13] = get_rsbox(Y[13]);
	Y[14] = get_rsbox(Y[14]);
	Y[15] = get_rsbox(Y[15]);
		
	/// AddRoundKey(Y,RK,0);
	Y[12] ^= RK[0][12];
	Y[13] ^= RK[0][13];
	Y[14] ^= RK[0][14];
	Y[15] ^= RK[0][15];
	
	/// Ultima
	return get_rsbox(Mult14(Y[12])^Mult11(Y[13])^Mult13(Y[14])^ Mult9(Y[15]));
}


void ConstructBiclique(byte C0[16],byte K0[16],byte Ci[256][16],byte Sj[256][16]){
	byte S0[16],*aux;
	int i,j;
	
	f_1(C0,K0,S0);
	
	for(i = 0; i < 256;i++){
		aux = XOR_matrix(K0,delta_k_i);
		f(S0,aux,Ci[i]);
		delta_k_i[8]++;
		delta_k_i[12]++;
		free(aux);
	}
	
	for(j = 0; j < 256;j++){
		aux = XOR_matrix(K0,nabla_k_j);
		f_1(C0,aux,Sj[j]);
		nabla_k_j[1]++;
		nabla_k_j[9]++;
		free(aux);
	}
	
	f_1(Ci[15],K0,S0);	
}

void GeneratePi(byte real_K[16],byte Pi[256][16],byte Ci[256][16]){
	int i;
	
	for(i = 0; i < 256;i++) AES_128_1(Ci[i],real_K,Pi[i]);
	
}

int MITM(byte K[16],byte Pi[256][16],byte Sj[256][16],byte real_K[16]){
	int i,j;
	byte testK[16],testS[16],RK[9][16];
	
	for(i = 0; i < 256; i++){
		for(j = 0; j < 256; j++){	
			memcpy(testK,K,16);
			testK[8] ^= i;
			testK[12]^= i;
			testK[1] ^= j;
			testK[9] ^= j;
			
			
			g(Pi[i],testK,testS);
			
			if(equals(testS,Sj[j])){
				BicliqueKeyExpansion_1(testK,RK);
				memcpy(real_K,RK[0],16);
				return 1;
			} 
			
		}
	}
	return 0;
}
	
void Precomputation(byte K[16],byte Pi[256][16],byte Sj[256][16]){
	int i,j;
	byte RK[9][16],real_RK[11][16],*aux;
	
	BicliqueKeyExpansion_1(K,real_RK);
	
	for(i = 0; i < 256;i++){		
		aux = XOR_matrix(K,delta_k_i);
		BicliqueKeyExpansion_1(aux,RKi[i]);
		delta_k_i[8]++;
		delta_k_i[12]++;
		vi[i][0] = r(Pi[i],fullforward[i],RKi[i]);
		free(aux);
	}	
	for(j = 0; j < 256;j++){
		aux = XOR_matrix(K,nabla_k_j);
		BicliqueKeyExpansion_1(aux,RKj[j]);
		nabla_k_j[1]++;
		nabla_k_j[9]++;
		t_1(Sj[j],fullbackward[j],RKj[j]);
		free(aux);
	}
}

int Recomputation(byte K[16],byte Pi[256][16],byte Sj[256][16],byte res[16]){
	int i,j,k,cont = 0,s=0;
	byte RK[5][16],RKaux[9][16],v;
	
	for(i = 0; i < 256; i++){
		for(j = 0; j < 256; j++){
			forward_keys(RKi[i],RK,j);
			vi[i][j] = forward_r(Pi[i],fullforward[i],RK);
		}
	}
	
	for(i = 0; i < 256; i++){
		for(j = 0; j < 256; j++){
			backward_keys(RKj[j],RK,i);
			v = backward_t_1(Sj[j],fullbackward[j],RK);
			
			if(is_candidate(v,i,j)){
				if(check_candidate(Pi[i],K,Sj[j],i,j)){
					memcpy(res,RKi[i][8],16);
					res[1] ^= j;
					res[9] ^= j;
					BicliqueKeyExpansion_1(res,RKaux);
					memcpy(res,RKaux[0],16);
					s = 1;
				}
			}
		}
	}
	
	return s;
}

int MP(byte K[16],byte Pi[256][16],byte Sj[256][16],byte real_K[16]){
	Precomputation(K,Pi,Sj);
	int res = Recomputation(K,Pi,Sj,real_K);
	return res;
}

int BicliqueAttack_MITM(byte C0[16],byte real_K[16], byte K0[16]){
	byte Ci[256][16],Sj[256][16],Pi[256][16],Kzin[16];
	
	ConstructBiclique(C0,K0,Ci,Sj);
	
	GeneratePi(real_K,Pi,Ci);
	
	int res = MITM(K0,Pi,Sj,Kzin);
	
	if(res) memcpy(K0,Kzin,16);
	
	return res;
}

int BicliqueAttack_MP(byte C0[16],byte real_K[16], byte K0[16]){
	byte b1,b2,S[16],Ci[256][16],Sj[256][16],Pi[256][16],Kzin[16],K1[16],K2[16],RK1[9][16],RK2[9][16],RK[5][16],pre[7],pre2[4][16];
	int res;
	
	ConstructBiclique(C0,K0,Ci,Sj);
	
	GeneratePi(real_K,Pi,Ci);
	
	res = MP(K0,Pi,Sj,Kzin);
	
	if(res) memcpy(K0,Kzin,16);
	
	return res;
}

