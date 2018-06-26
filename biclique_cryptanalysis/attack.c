#include "aes.h"
#include "biclique.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "string.h"
#include "math.h"

static int LFSR_key[4] = {0x9AF5AB20,0x356A9CAB,0x0350F6CF,0x8C3BBE48};
byte lfsr[16];
byte lfsr_helper[16];
byte generator_helper[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
byte generator_helper2[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

byte K[16] = {0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c};

int calculate_percent(double percent,int i){
	double value = ((double)i)/ percent;
	int res = floor(value);
	if(value-res > 0.5) return res++;
	return res;
}

void fix_bytes(byte K2[16],int f){///assings the first f positions of K2 to be equal to the first f positions of global K
	if(f < 0 || f > 16) return;
	int i;	
	memcpy(K2,K,f);
	for(i = f; i < 16; i++) K2[i] = 0;
}

void biclique_fix_bytes(byte K[16],byte K2[16],int f){///assings the first f positions of K2 to be equal to the first f positions of K 
													  ///with exception of positions 1, 8, 9 and 12
	if(f < 0 || f > 14) return;
	int i,z = 16-f-2;	
	
	memcpy(K2,K,16);
	K2[8] ^= K2[12];
	K2[9] ^= K2[1];

	K2[1] = 0;
	K2[12] = 0;	
	
	if(z > 0){
		z--;
		K2[9] = 0;
		
		if(z > 0){
			z--;
			K2[8] = 0;
		}
	}
	
	for(i = 15; i >= 16-z; i--){
		if(i != 1 && i != 8 && i != 9 && i != 12) K2[i] = 0;
		else z++;
	} 
}

void biclique_copy_bytes(byte *K, byte *K2, int f){///assings the first 'f' positions of 'K' (starting
													///from 'ini') to be equal to the first 'f' positions of 'K2'
													///excluding positions 1 and 12
	if(f <= 0 || f > 14) return;
	
	int i,j = 0,ini,aux;	
	
	if(f > 0){
		f--;
		K[9] = K2[j++];
		if(f > 0){
			f--;
			K[8] = K2[j++];
		}
	}

	
	ini = 16-f;
	
	for(i = 15; i >=ini; i--){	
		if(i != 1 && i != 8 && i != 9 && i != 12) K[i] = K2[j++];
		else ini--;
	}
}

void int_to_char(int *vet1, byte *vet2,int size,int num){/// turns 'num' bytes of int 'vet1' of size 'size' into 'num' 
														/// bytes of char 'vet2' of size 'size'*4
	if(size < 0) return;
	if(num < 1 || num > 16) return;
	
	int i,j = -1;
	
	for(i = 0; i < num; i++){
		if(j == -1) j = 1;
		else j = -1;
		vet2[i] = (vet1[i/4]>>((i+(j+2)%4)*8))&0xff;
	}
}

void generate_1_bytes(byte *lfsr,int j){/// generates 1 byte from '*lfsr'

	byte aux,aux1,temp,vet[4] = {0,0,0,0};
	int i;
	if(generator_helper2[j] == 0xfe){
		*lfsr = 0xff;
		generator_helper2[j]++;
		return;
	}
	
	if(generator_helper2[j] == 0xff){
		*lfsr = lfsr_helper[j];
		generator_helper2[j]++;
		return;
	}
	
	for(i = 0; i < 256; i++){
		
		aux = *lfsr & 1;
		*lfsr = (*lfsr>>1)&0x7f;
		*lfsr = *lfsr^(aux<<7);
		aux = (aux+1)%2;
		*lfsr = *lfsr^(aux<<5);
		*lfsr = *lfsr^(aux<<4);
		*lfsr = *lfsr^(aux<<3);
	}
	generator_helper2[j]++;
}

void generate_16_bytes(byte X[16]){/// generates a random state
	int i,j, vet[4] = {0,0,0,0}, temp[3] = {0,0,0}, aux;
	
	int_to_char(LFSR_key,X,4,16);
	
    for(i = 0; i < 128; i++){
 
        temp[0] = (LFSR_key[0] & 1)<<31;
        temp[1] = (LFSR_key[1] & 1)<<31;
        temp[2] = (LFSR_key[2] & 1)<<31;
 
        ///D^256 + D^254 + D^251 + D^246 + 1
        vet[0] =  LFSR_key[0]&0x80000000;
        vet[1] = (LFSR_key[0]<<2)&0x80000000;
        vet[2] = (LFSR_key[0]<<5)&0x80000000;
        vet[3] = (LFSR_key[0]<<10)&0x80000000;
		
        aux = vet[0]^vet[1]^vet[2]^vet[3];
 
        for(j = 0; j < 4; j++){
            if (j == 0)
                LFSR_key[j] = ((LFSR_key[j]>>1)&0x7fffffff) ^ aux;
            else
                LFSR_key[j] = ((LFSR_key[j]>>1)&0x7fffffff) ^ temp[j-1];
        }
    }
}

void Generate_random_bytes(byte *X, int num){/// generates 'num' bytes from 'lfsr' and puts them in 'X'
	
	if(num < 1 || num > 16) return;
	
	int i;
    
	memcpy(X,lfsr,num);
	
	for(i = 0; i < num; i++){
		generate_1_bytes(&lfsr[i],i);
		generator_helper[i]++;
		if(generator_helper[i] != 0) break;
	}
}
 
void Biclique_MITM(int d, int f){/// K key, d dimension of biclique, f # fixed bytes
	byte K0[16],P[16],C[16],C0[16],*aux,RK[11][16];
	
	aux = (byte *)malloc(sizeof(byte)*(16-f));
	KeyExpansion(K,RK);
	biclique_fix_bytes(RK[8],K0,f);

	generate_16_bytes(P);
	generate_16_bytes(C0);
	
	AES_128(P,K,C);
	
	int k = 8*(16-f),j,p = 128,s = 0	;
	unsigned long long i, total = (unsigned long long)pow(2,k-(2*d));
	double percent = ((double)total)/100;

	if(total == 1){
		total = 128;
		p = 1;
	}
	
	printf("Biclique_MITM - #grupos de chaves: %llx\n",total);
	
	int start = clock();	 
	for(j = 0; j < p; j++){
		printf("\r%.0f%%",percent);
		for(i = 0; i < total/128; i++){
			if(BicliqueAttack_MITM(C0,K,K0)){
				printf("\rSuccess\n");
				printState(K0);
				s=1;
			}
			
			Generate_random_bytes(aux,14-f);
			
			biclique_copy_bytes(K0,aux,14-f);
		}
		percent += 0.78125;
	}
	printf("\r100%%");
	int end = clock();	 
	printf("\ntime: %f\n",((double)(end-start))/1000);
	FILE *file = fopen("result-MITM.txt","a");
	if(s) fprintf(file,"success\n");
	fprintf(file,"v1.0:%llx,	%d\n",i,f);
	fprintf(file,"time: %fs\n",((double)(end-start))/1000);
	fprintf(file,"time: %fm\n",((double)(end-start))/1000/60);
	fprintf(file,"time: %fh\n",((double)(end-start))/1000/3600);
	fclose(file);
	free(aux);
}

void Biclique_MP(int d, int f){/// K key, d dimension of biclique, f # fixed bytes
	byte K0[16],P[16],C[16],C0[16],*aux,RK[11][16];
	
	aux = (byte *)malloc(sizeof(byte)*(16-f));
	KeyExpansion(K,RK);
	biclique_fix_bytes(RK[8],K0,f);

	generate_16_bytes(P);
	generate_16_bytes(C0);
	
	AES_128(P,K,C);
	
	int k = 8*(16-f),s=0,j,p = 128;	
	unsigned long long i, total = (unsigned long long)pow(2,k-(2*d));
	double percent = 0;
	
	if(total == 1){
		total = 128;
		p = 1;
	}
	
	printf("Biclique_MP - #grupos de chaves: %llx\n",total);

	int start = clock();	
	for(j = 0; j < p; j++){
		printf("\r%.0f%%",percent);
		for(i = 0; i < total/128; i++){
			if(BicliqueAttack_MP(C0,K,K0)){
				printf("\rSuccess\n");
				printState(K0);
				s=1;
			}
			
			Generate_random_bytes(aux,14-f);
			
			biclique_copy_bytes(K0,aux,14-f);
		}
		percent += 0.78125;
	}
	printf("\r100%%");
	int end = clock();	
	printf("\ntime: %f\n",((double)(end-start))/1000);
	FILE *file = fopen("result-MP.txt","a");
	fprintf(file,"v1.0:%llx,	%d\n",i,f);
	if(s) fprintf(file,"success\n");
	fprintf(file,"time: %fs\n",((double)(end-start))/1000);
	fprintf(file,"time: %fm\n",((double)(end-start))/1000/60);
	fprintf(file,"time: %fh\n",((double)(end-start))/1000/3600);
	fclose(file);
	free(aux);
}

void Brute_force(int f){/// K key, f # fixed bytes
	byte K2[16],P[16],C[16],C2[16],*aux;

	aux = (byte *)malloc(sizeof(byte)*(16-f));
	fix_bytes(K2,f);
	generate_16_bytes(P);
	
	AES_128(P,K,C);
	
	int k = 8*(16-f),s = 0,j,p = 128;
	unsigned long long i, total = (unsigned long long)pow(2,k);
	double percent = 0;
	
	if(total == 1){
		total = 128;
		p = 1;
	}
	
	printf("Brute Force - #grupos de chaves: %llx\n",total);
	
	int start = clock();
	for(j = 0; j < p; j++){	
		printf("\r%.0f%%",percent);
		for(i = 0; i < total/128; i++){
			
			AES_128(P,K2,C2);
			
			if(equals(C,C2)){
				printf("\rSuccess!!\n");
				printState(K2);

				s = i;
			}
			
			Generate_random_bytes(aux,16-f);
			memcpy(K2+f,aux,(16-f));
		}
		percent += 0.78125;
	}
	
	printf("\r100%%");
	int end = clock();	 
	printf("    %f\n",((double)(end-start))/1000);
	FILE *file = fopen("result_BFs.txt","a");
	fprintf(file,"v1.0:%llx,	%d\n",i,f);
	if(s) fprintf(file,"noice\n");
	fprintf(file,"time: %fs\n",((double)(end-start))/1000);
	fprintf(file,"time: %fm\n",((double)(end-start))/1000/60);
	fprintf(file,"time: %fh\n",((double)(end-start))/1000/3600);
	fclose(file);
	free(aux);
}

int main(){
	int_to_char(LFSR_key,lfsr,4,16);
	memcpy(lfsr_helper,lfsr,16);
	int i,j,seed = 3;
	byte K2[16] = {0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0};
	byte RK[11][16];
	
	system("cls");	
	
	for(i = 0; i++; i < seed) generate_16_bytes(K2);
	
	memcpy(K,K2,16);
	printf("K\n");
	printState(K);
	KeyExpansion(K,RK);
	printState(RK[8]);

	// Biclique_MITM(8,13);
	Biclique_MP(8,14);
	// Brute_force(13);
	
	return 0;
}

