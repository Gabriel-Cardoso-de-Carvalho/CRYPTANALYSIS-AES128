#include "aes.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"

/// The lookup-tables are marked const so they can be placed in read-only storage instead of RAM
/// The numbers below can be computed dynamically trading ROM for RAM - 
/// This can be useful in (embedded) bootloader applications, where ROM is often limited.
static const byte sbox[256] = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

static const byte rsbox[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

/// The round constant word array, Rcon[i], contains the values given by 
/// x to the power (i-1) being powers of x (x is denoted as {02}) in the field GF(2^8)
static const byte Rcon[10] = {
  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36 };

static const byte M2[256]={
	0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x10, 0x12, 0x14, 0x16, 0x18, 0x1a, 0x1c, 0x1e, 0x20, 
	0x22, 0x24, 0x26, 0x28, 0x2a, 0x2c, 0x2e, 0x30, 0x32, 0x34, 0x36, 0x38, 0x3a, 0x3c, 0x3e, 0x40, 
	0x42, 0x44, 0x46, 0x48, 0x4a, 0x4c, 0x4e, 0x50, 0x52, 0x54, 0x56, 0x58, 0x5a, 0x5c, 0x5e, 0x60, 
	0x62, 0x64, 0x66, 0x68, 0x6a, 0x6c, 0x6e, 0x70, 0x72, 0x74, 0x76, 0x78, 0x7a, 0x7c, 0x7e, 0x80, 
	0x82, 0x84, 0x86, 0x88, 0x8a, 0x8c, 0x8e, 0x90, 0x92, 0x94, 0x96, 0x98, 0x9a, 0x9c, 0x9e, 0xa0, 
	0xa2, 0xa4, 0xa6, 0xa8, 0xaa, 0xac, 0xae, 0xb0, 0xb2, 0xb4, 0xb6, 0xb8, 0xba, 0xbc, 0xbe, 0xc0, 
	0xc2, 0xc4, 0xc6, 0xc8, 0xca, 0xcc, 0xce, 0xd0, 0xd2, 0xd4, 0xd6, 0xd8, 0xda, 0xdc, 0xde, 0xe0, 
	0xe2, 0xe4, 0xe6, 0xe8, 0xea, 0xec, 0xee, 0xf0, 0xf2, 0xf4, 0xf6, 0xf8, 0xfa, 0xfc, 0xfe, 0x1b, 
	0x19, 0x1f, 0x1d, 0x13, 0x11, 0x17, 0x15, 0xb, 0x9, 0xf, 0xd, 0x3, 0x1, 0x7, 0x5, 0x3b, 
	0x39, 0x3f, 0x3d, 0x33, 0x31, 0x37, 0x35, 0x2b, 0x29, 0x2f, 0x2d, 0x23, 0x21, 0x27, 0x25, 0x5b, 
	0x59, 0x5f, 0x5d, 0x53, 0x51, 0x57, 0x55, 0x4b, 0x49, 0x4f, 0x4d, 0x43, 0x41, 0x47, 0x45, 0x7b, 
	0x79, 0x7f, 0x7d, 0x73, 0x71, 0x77, 0x75, 0x6b, 0x69, 0x6f, 0x6d, 0x63, 0x61, 0x67, 0x65, 0x9b, 
	0x99, 0x9f, 0x9d, 0x93, 0x91, 0x97, 0x95, 0x8b, 0x89, 0x8f, 0x8d, 0x83, 0x81, 0x87, 0x85, 0xbb, 
	0xb9, 0xbf, 0xbd, 0xb3, 0xb1, 0xb7, 0xb5, 0xab, 0xa9, 0xaf, 0xad, 0xa3, 0xa1, 0xa7, 0xa5, 0xdb, 
	0xd9, 0xdf, 0xdd, 0xd3, 0xd1, 0xd7, 0xd5, 0xcb, 0xc9, 0xcf, 0xcd, 0xc3, 0xc1, 0xc7, 0xc5, 0xfb, 
	0xf9, 0xff, 0xfd, 0xf3, 0xf1, 0xf7, 0xf5, 0xeb, 0xe9, 0xef, 0xed, 0xe3, 0xe1, 0xe7, 0xe5 
};
  
static const byte M3[256]={
	0x0, 0x3, 0x6, 0x5, 0xc, 0xf, 0xa, 0x9, 0x18, 0x1b, 0x1e, 0x1d, 0x14, 0x17, 0x12, 0x11, 0x30, 
	0x33, 0x36, 0x35, 0x3c, 0x3f, 0x3a, 0x39, 0x28, 0x2b, 0x2e, 0x2d, 0x24, 0x27, 0x22, 0x21, 0x60, 
	0x63, 0x66, 0x65, 0x6c, 0x6f, 0x6a, 0x69, 0x78, 0x7b, 0x7e, 0x7d, 0x74, 0x77, 0x72, 0x71, 0x50, 
	0x53, 0x56, 0x55, 0x5c, 0x5f, 0x5a, 0x59, 0x48, 0x4b, 0x4e, 0x4d, 0x44, 0x47, 0x42, 0x41, 0xc0, 
	0xc3, 0xc6, 0xc5, 0xcc, 0xcf, 0xca, 0xc9, 0xd8, 0xdb, 0xde, 0xdd, 0xd4, 0xd7, 0xd2, 0xd1, 0xf0, 
	0xf3, 0xf6, 0xf5, 0xfc, 0xff, 0xfa, 0xf9, 0xe8, 0xeb, 0xee, 0xed, 0xe4, 0xe7, 0xe2, 0xe1, 0xa0, 
	0xa3, 0xa6, 0xa5, 0xac, 0xaf, 0xaa, 0xa9, 0xb8, 0xbb, 0xbe, 0xbd, 0xb4, 0xb7, 0xb2, 0xb1, 0x90, 
	0x93, 0x96, 0x95, 0x9c, 0x9f, 0x9a, 0x99, 0x88, 0x8b, 0x8e, 0x8d, 0x84, 0x87, 0x82, 0x81, 0x9b, 
	0x98, 0x9d, 0x9e, 0x97, 0x94, 0x91, 0x92, 0x83, 0x80, 0x85, 0x86, 0x8f, 0x8c, 0x89, 0x8a, 0xab, 
	0xa8, 0xad, 0xae, 0xa7, 0xa4, 0xa1, 0xa2, 0xb3, 0xb0, 0xb5, 0xb6, 0xbf, 0xbc, 0xb9, 0xba, 0xfb, 
	0xf8, 0xfd, 0xfe, 0xf7, 0xf4, 0xf1, 0xf2, 0xe3, 0xe0, 0xe5, 0xe6, 0xef, 0xec, 0xe9, 0xea, 0xcb, 
	0xc8, 0xcd, 0xce, 0xc7, 0xc4, 0xc1, 0xc2, 0xd3, 0xd0, 0xd5, 0xd6, 0xdf, 0xdc, 0xd9, 0xda, 0x5b, 
	0x58, 0x5d, 0x5e, 0x57, 0x54, 0x51, 0x52, 0x43, 0x40, 0x45, 0x46, 0x4f, 0x4c, 0x49, 0x4a, 0x6b, 
	0x68, 0x6d, 0x6e, 0x67, 0x64, 0x61, 0x62, 0x73, 0x70, 0x75, 0x76, 0x7f, 0x7c, 0x79, 0x7a, 0x3b, 
	0x38, 0x3d, 0x3e, 0x37, 0x34, 0x31, 0x32, 0x23, 0x20, 0x25, 0x26, 0x2f, 0x2c, 0x29, 0x2a, 0xb, 
	0x8, 0xd, 0xe, 0x7, 0x4, 0x1, 0x2, 0x13, 0x10, 0x15, 0x16, 0x1f, 0x1c, 0x19, 0x1a
};

static const byte M9[256]={
	0x0, 0x9, 0x12, 0x1b, 0x24, 0x2d, 0x36, 0x3f, 0x48, 0x41, 0x5a, 0x53, 0x6c, 0x65, 0x7e, 0x77, 0x90, 
	0x99, 0x82, 0x8b, 0xb4, 0xbd, 0xa6, 0xaf, 0xd8, 0xd1, 0xca, 0xc3, 0xfc, 0xf5, 0xee, 0xe7, 0x3b, 
	0x32, 0x29, 0x20, 0x1f, 0x16, 0xd, 0x4, 0x73, 0x7a, 0x61, 0x68, 0x57, 0x5e, 0x45, 0x4c, 0xab, 
	0xa2, 0xb9, 0xb0, 0x8f, 0x86, 0x9d, 0x94, 0xe3, 0xea, 0xf1, 0xf8, 0xc7, 0xce, 0xd5, 0xdc, 0x76, 
	0x7f, 0x64, 0x6d, 0x52, 0x5b, 0x40, 0x49, 0x3e, 0x37, 0x2c, 0x25, 0x1a, 0x13, 0x8, 0x1, 0xe6, 
	0xef, 0xf4, 0xfd, 0xc2, 0xcb, 0xd0, 0xd9, 0xae, 0xa7, 0xbc, 0xb5, 0x8a, 0x83, 0x98, 0x91, 0x4d, 
	0x44, 0x5f, 0x56, 0x69, 0x60, 0x7b, 0x72, 0x5, 0xc, 0x17, 0x1e, 0x21, 0x28, 0x33, 0x3a, 0xdd, 
	0xd4, 0xcf, 0xc6, 0xf9, 0xf0, 0xeb, 0xe2, 0x95, 0x9c, 0x87, 0x8e, 0xb1, 0xb8, 0xa3, 0xaa, 0xec, 
	0xe5, 0xfe, 0xf7, 0xc8, 0xc1, 0xda, 0xd3, 0xa4, 0xad, 0xb6, 0xbf, 0x80, 0x89, 0x92, 0x9b, 0x7c, 
	0x75, 0x6e, 0x67, 0x58, 0x51, 0x4a, 0x43, 0x34, 0x3d, 0x26, 0x2f, 0x10, 0x19, 0x2, 0xb, 0xd7, 
	0xde, 0xc5, 0xcc, 0xf3, 0xfa, 0xe1, 0xe8, 0x9f, 0x96, 0x8d, 0x84, 0xbb, 0xb2, 0xa9, 0xa0, 0x47, 
	0x4e, 0x55, 0x5c, 0x63, 0x6a, 0x71, 0x78, 0xf, 0x6, 0x1d, 0x14, 0x2b, 0x22, 0x39, 0x30, 0x9a, 
	0x93, 0x88, 0x81, 0xbe, 0xb7, 0xac, 0xa5, 0xd2, 0xdb, 0xc0, 0xc9, 0xf6, 0xff, 0xe4, 0xed, 0xa, 
	0x3, 0x18, 0x11, 0x2e, 0x27, 0x3c, 0x35, 0x42, 0x4b, 0x50, 0x59, 0x66, 0x6f, 0x74, 0x7d, 0xa1, 
	0xa8, 0xb3, 0xba, 0x85, 0x8c, 0x97, 0x9e, 0xe9, 0xe0, 0xfb, 0xf2, 0xcd, 0xc4, 0xdf, 0xd6, 0x31, 
	0x38, 0x23, 0x2a, 0x15, 0x1c, 0x7, 0xe, 0x79, 0x70, 0x6b, 0x62, 0x5d, 0x54, 0x4f, 0x46 
};

static const byte M11[256]={
	0x0, 0xb, 0x16, 0x1d, 0x2c, 0x27, 0x3a, 0x31, 0x58, 0x53, 0x4e, 0x45, 0x74, 0x7f, 0x62, 0x69, 0xb0, 
	0xbb, 0xa6, 0xad, 0x9c, 0x97, 0x8a, 0x81, 0xe8, 0xe3, 0xfe, 0xf5, 0xc4, 0xcf, 0xd2, 0xd9, 0x7b, 
	0x70, 0x6d, 0x66, 0x57, 0x5c, 0x41, 0x4a, 0x23, 0x28, 0x35, 0x3e, 0xf, 0x4, 0x19, 0x12, 0xcb, 
	0xc0, 0xdd, 0xd6, 0xe7, 0xec, 0xf1, 0xfa, 0x93, 0x98, 0x85, 0x8e, 0xbf, 0xb4, 0xa9, 0xa2, 0xf6, 
	0xfd, 0xe0, 0xeb, 0xda, 0xd1, 0xcc, 0xc7, 0xae, 0xa5, 0xb8, 0xb3, 0x82, 0x89, 0x94, 0x9f, 0x46, 
	0x4d, 0x50, 0x5b, 0x6a, 0x61, 0x7c, 0x77, 0x1e, 0x15, 0x8, 0x3, 0x32, 0x39, 0x24, 0x2f, 0x8d, 
	0x86, 0x9b, 0x90, 0xa1, 0xaa, 0xb7, 0xbc, 0xd5, 0xde, 0xc3, 0xc8, 0xf9, 0xf2, 0xef, 0xe4, 0x3d, 
	0x36, 0x2b, 0x20, 0x11, 0x1a, 0x7, 0xc, 0x65, 0x6e, 0x73, 0x78, 0x49, 0x42, 0x5f, 0x54, 0xf7, 
	0xfc, 0xe1, 0xea, 0xdb, 0xd0, 0xcd, 0xc6, 0xaf, 0xa4, 0xb9, 0xb2, 0x83, 0x88, 0x95, 0x9e, 0x47, 
	0x4c, 0x51, 0x5a, 0x6b, 0x60, 0x7d, 0x76, 0x1f, 0x14, 0x9, 0x2, 0x33, 0x38, 0x25, 0x2e, 0x8c, 
	0x87, 0x9a, 0x91, 0xa0, 0xab, 0xb6, 0xbd, 0xd4, 0xdf, 0xc2, 0xc9, 0xf8, 0xf3, 0xee, 0xe5, 0x3c, 
	0x37, 0x2a, 0x21, 0x10, 0x1b, 0x6, 0xd, 0x64, 0x6f, 0x72, 0x79, 0x48, 0x43, 0x5e, 0x55, 0x1, 
	0xa, 0x17, 0x1c, 0x2d, 0x26, 0x3b, 0x30, 0x59, 0x52, 0x4f, 0x44, 0x75, 0x7e, 0x63, 0x68, 0xb1, 
	0xba, 0xa7, 0xac, 0x9d, 0x96, 0x8b, 0x80, 0xe9, 0xe2, 0xff, 0xf4, 0xc5, 0xce, 0xd3, 0xd8, 0x7a, 
	0x71, 0x6c, 0x67, 0x56, 0x5d, 0x40, 0x4b, 0x22, 0x29, 0x34, 0x3f, 0xe, 0x5, 0x18, 0x13, 0xca, 
	0xc1, 0xdc, 0xd7, 0xe6, 0xed, 0xf0, 0xfb, 0x92, 0x99, 0x84, 0x8f, 0xbe, 0xb5, 0xa8, 0xa3,
};

static const byte M13[256]={
	0x0, 0xd, 0x1a, 0x17, 0x34, 0x39, 0x2e, 0x23, 0x68, 0x65, 0x72, 0x7f, 0x5c, 0x51, 0x46, 0x4b, 0xd0, 
	0xdd, 0xca, 0xc7, 0xe4, 0xe9, 0xfe, 0xf3, 0xb8, 0xb5, 0xa2, 0xaf, 0x8c, 0x81, 0x96, 0x9b, 0xbb, 
	0xb6, 0xa1, 0xac, 0x8f, 0x82, 0x95, 0x98, 0xd3, 0xde, 0xc9, 0xc4, 0xe7, 0xea, 0xfd, 0xf0, 0x6b, 
	0x66, 0x71, 0x7c, 0x5f, 0x52, 0x45, 0x48, 0x3, 0xe, 0x19, 0x14, 0x37, 0x3a, 0x2d, 0x20, 0x6d, 
	0x60, 0x77, 0x7a, 0x59, 0x54, 0x43, 0x4e, 0x5, 0x8, 0x1f, 0x12, 0x31, 0x3c, 0x2b, 0x26, 0xbd, 
	0xb0, 0xa7, 0xaa, 0x89, 0x84, 0x93, 0x9e, 0xd5, 0xd8, 0xcf, 0xc2, 0xe1, 0xec, 0xfb, 0xf6, 0xd6, 
	0xdb, 0xcc, 0xc1, 0xe2, 0xef, 0xf8, 0xf5, 0xbe, 0xb3, 0xa4, 0xa9, 0x8a, 0x87, 0x90, 0x9d, 0x6, 
	0xb, 0x1c, 0x11, 0x32, 0x3f, 0x28, 0x25, 0x6e, 0x63, 0x74, 0x79, 0x5a, 0x57, 0x40, 0x4d, 0xda, 
	0xd7, 0xc0, 0xcd, 0xee, 0xe3, 0xf4, 0xf9, 0xb2, 0xbf, 0xa8, 0xa5, 0x86, 0x8b, 0x9c, 0x91, 0xa, 
	0x7, 0x10, 0x1d, 0x3e, 0x33, 0x24, 0x29, 0x62, 0x6f, 0x78, 0x75, 0x56, 0x5b, 0x4c, 0x41, 0x61, 
	0x6c, 0x7b, 0x76, 0x55, 0x58, 0x4f, 0x42, 0x9, 0x4, 0x13, 0x1e, 0x3d, 0x30, 0x27, 0x2a, 0xb1, 
	0xbc, 0xab, 0xa6, 0x85, 0x88, 0x9f, 0x92, 0xd9, 0xd4, 0xc3, 0xce, 0xed, 0xe0, 0xf7, 0xfa, 0xb7, 
	0xba, 0xad, 0xa0, 0x83, 0x8e, 0x99, 0x94, 0xdf, 0xd2, 0xc5, 0xc8, 0xeb, 0xe6, 0xf1, 0xfc, 0x67, 
	0x6a, 0x7d, 0x70, 0x53, 0x5e, 0x49, 0x44, 0xf, 0x2, 0x15, 0x18, 0x3b, 0x36, 0x21, 0x2c, 0xc, 
	0x1, 0x16, 0x1b, 0x38, 0x35, 0x22, 0x2f, 0x64, 0x69, 0x7e, 0x73, 0x50, 0x5d, 0x4a, 0x47, 0xdc, 
	0xd1, 0xc6, 0xcb, 0xe8, 0xe5, 0xf2, 0xff, 0xb4, 0xb9, 0xae, 0xa3, 0x80, 0x8d, 0x9a, 0x97
};

static const byte M14[256]={
	0x0, 0xe, 0x1c, 0x12, 0x38, 0x36, 0x24, 0x2a, 0x70, 0x7e, 0x6c, 0x62, 0x48, 0x46, 0x54, 0x5a, 0xe0, 
	0xee, 0xfc, 0xf2, 0xd8, 0xd6, 0xc4, 0xca, 0x90, 0x9e, 0x8c, 0x82, 0xa8, 0xa6, 0xb4, 0xba, 0xdb, 
	0xd5, 0xc7, 0xc9, 0xe3, 0xed, 0xff, 0xf1, 0xab, 0xa5, 0xb7, 0xb9, 0x93, 0x9d, 0x8f, 0x81, 0x3b, 
	0x35, 0x27, 0x29, 0x3, 0xd, 0x1f, 0x11, 0x4b, 0x45, 0x57, 0x59, 0x73, 0x7d, 0x6f, 0x61, 0xad, 
	0xa3, 0xb1, 0xbf, 0x95, 0x9b, 0x89, 0x87, 0xdd, 0xd3, 0xc1, 0xcf, 0xe5, 0xeb, 0xf9, 0xf7, 0x4d, 
	0x43, 0x51, 0x5f, 0x75, 0x7b, 0x69, 0x67, 0x3d, 0x33, 0x21, 0x2f, 0x5, 0xb, 0x19, 0x17, 0x76, 
	0x78, 0x6a, 0x64, 0x4e, 0x40, 0x52, 0x5c, 0x6, 0x8, 0x1a, 0x14, 0x3e, 0x30, 0x22, 0x2c, 0x96, 
	0x98, 0x8a, 0x84, 0xae, 0xa0, 0xb2, 0xbc, 0xe6, 0xe8, 0xfa, 0xf4, 0xde, 0xd0, 0xc2, 0xcc, 0x41, 
	0x4f, 0x5d, 0x53, 0x79, 0x77, 0x65, 0x6b, 0x31, 0x3f, 0x2d, 0x23, 0x9, 0x7, 0x15, 0x1b, 0xa1, 
	0xaf, 0xbd, 0xb3, 0x99, 0x97, 0x85, 0x8b, 0xd1, 0xdf, 0xcd, 0xc3, 0xe9, 0xe7, 0xf5, 0xfb, 0x9a, 
	0x94, 0x86, 0x88, 0xa2, 0xac, 0xbe, 0xb0, 0xea, 0xe4, 0xf6, 0xf8, 0xd2, 0xdc, 0xce, 0xc0, 0x7a, 
	0x74, 0x66, 0x68, 0x42, 0x4c, 0x5e, 0x50, 0xa, 0x4, 0x16, 0x18, 0x32, 0x3c, 0x2e, 0x20, 0xec, 
	0xe2, 0xf0, 0xfe, 0xd4, 0xda, 0xc8, 0xc6, 0x9c, 0x92, 0x80, 0x8e, 0xa4, 0xaa, 0xb8, 0xb6, 0xc, 
	0x2, 0x10, 0x1e, 0x34, 0x3a, 0x28, 0x26, 0x7c, 0x72, 0x60, 0x6e, 0x44, 0x4a, 0x58, 0x56, 0x37, 
	0x39, 0x2b, 0x25, 0xf, 0x1, 0x13, 0x1d, 0x47, 0x49, 0x5b, 0x55, 0x7f, 0x71, 0x63, 0x6d, 0xd7, 
	0xd9, 0xcb, 0xc5, 0xef, 0xe1, 0xf3, 0xfd, 0xa7, 0xa9, 0xbb, 0xb5, 0x9f, 0x91, 0x83, 0x8d
};
  
byte *XOR_matrix(byte M1[16],byte M2[16]){
	byte *M = (byte *)malloc(16*sizeof(byte));
	
	int i;
	for(i = 0; i < 16;i++){
		M[i] = M1[i]^M2[i];
	}
	return M;
}
  
byte Mult2(byte a){
	/// return ((a>>7)&1)*0x1B ^ (a<<1);
	return M2[a];
}

byte Mult3(byte a){
	/// return Mult2(a)^a;
	return M3[a];
}

byte Mult9(byte a){
	/// return Mult2(Mult2(Mult2(a)))^a;
	return M9[a];
}

byte Mult11(byte a){
	/// return Mult2(Mult2(Mult2(a)))^Mult2(a)^a;
	return M11[a];
}

byte Mult13(byte a){
	/// return Mult2(Mult2(Mult2(a)))^Mult2(Mult2(a))^a;
	return M13[a];
}

byte Mult14(byte a){
	/// return Mult2(Mult2(Mult2(a)))^Mult2(Mult2(a))^Mult2(a);
	return M14[a];
}
  
byte get_sbox(byte b){
	return sbox[b];
}  

byte get_rsbox(byte b){
	return rsbox[b];
}  
  
void printState(byte X[16]){
	int i,j;
	for(i = 0; i < 4;i++){
		for(j = 0; j < 4;j++) printf("%x ",(X[i+j*4])&0xff);
		printf("\n");
	}
	printf("\n");
}  

void print2States(byte X[16], byte Y[16]){
	int i,j;
	for(i = 0; i < 4;i++){
		for(j = 0; j < 4;j++){
			printf("%x ",(X[i+j*4])&0xff);
		}
		printf("        ");
		for(j = 0; j < 4;j++){
			printf("%x ",(Y[i+j*4])&0xff);
		}
		printf("\n");
	}
	printf("\n");
}

int equals(byte C[16],byte C2[16]){
	int i;
	
	for(i = 0; i < 16; i++){
		if(C[i] != C2[i]){
			return 0;
		}
	}
	return 1;
}

void KeySchedule(byte K[16],int round){
	if(round >=10) return;
	byte aux[4];
	int i;
	
	for(i = 0; i < 3;i++) aux[i] = sbox[K[i+13]];
	aux[3] = sbox[K[12]];
	
	K[0] = K[0] ^ aux[0] ^ Rcon[round];
	for(i = 1; i < 4;i++) K[i] = K[i] ^ aux[i];
	
	for(i = 4; i < 16;i++){
		K[i] = K[i]^K[i-4];
	}
}

void KeySchedule_1(byte K[16],int round){
	if(round >=10) return;
	byte aux[4];
	int i;
	
	for(i = 15; i >= 4;i--){
		K[i] = K[i]^K[i-4];
	}
		
	for(i = 0; i < 3;i++) aux[i] = sbox[K[i+13]];
	aux[3] = sbox[K[12]];
	
	K[0] = K[0] ^ aux[0] ^ Rcon[round];
	for(i = 1; i < 4;i++) K[i] = K[i] ^ aux[i];
	

}

void KeyExpansion(byte K[16],byte RK[11][16]){
	int i;
	
	
	for(i = 0; i < 10;i++){
		memcpy(RK[i],K,16);
		KeySchedule(K,i);
	}
	memcpy(RK[10],K,16);
	memcpy(K,RK[0],16);

}

void AddRoundKey(byte X[16], byte RK[11][16],int round){
	int i;
	if(round >=0 && round < 11){
		for(i = 0; i < 16;i++) X[i] = X[i]^RK[round][i];
	}
	
}

void AddRoundKey_1(byte X[16], byte K[16],int round){
	int i;
	
	for(i = 0; i < 16;i++) X[i] = X[i]^K[i];
	if(round >=0 && round < 10)
		KeySchedule(K,round);
	
}  
  
void SubBytes(byte X[16] ){
	int i;
	for(i = 0; i < 16;i++) X[i] = sbox[X[i]];
}
 
void SubBytes_1(byte X[16] ){
	int i;
	for(i = 0; i < 16;i++) X[i] = rsbox[X[i]];
} 
  
void ShiftRows(byte X[16] ){
	int i,j,temp[4];
	
	for(i = 0; i < 4;i++){
		
		temp[0] = X[i];
		temp[1] = X[4+i];
		temp[2] = X[8+i];
		temp[3] = X[12+i];
		
		for(j = 0; j < 4;j++){
			X[(j*4+i+(12*i))%16] = temp[j];
		}
	}
	
}

void ShiftRows_1(byte X[16] ){
	int i,j,temp[4],aux;
	
	for(i = 0; i < 4;i++){
		
		temp[0] = X[i];
		temp[1] = X[4+i];
		temp[2] = X[8+i];
		temp[3] = X[12+i];
		
		for(j = 0; j < 4;j++){
			aux = (j*4+i-(12*i))%16;
			if(aux < 0) aux += 16;
			X[aux] = temp[j];
		}
	}
	
}
 
void MixColumns(byte X[16]){
	int i,d[4];
	
	for(i = 0; i < 4;i++){
		d[0] = Mult2(X[i*4])^Mult3(X[i*4+1])^X[i*4+2]		^X[i*4+3];
		d[1] = X[i*4]	    ^Mult2(X[i*4+1])^Mult3(X[i*4+2])^X[i*4+3];
		d[2] = X[i*4]		^X[i*4+1]		^Mult2(X[i*4+2])^Mult3(X[i*4+3]);
		d[3] = Mult3(X[i*4])^X[i*4+1]		^X[i*4+2]		^Mult2(X[i*4+3]);
		
		X[i*4]   = d[0];
		X[i*4+1] = d[1];
		X[i*4+2] = d[2];
		X[i*4+3] = d[3];
	}
}

void MixColumns_1(byte X[16]){
	int i,d[4];
	
	for(i = 0; i < 4;i++){
		d[0] = Mult14(X[i*4])^Mult11(X[i*4+1])^Mult13(X[i*4+2])^ Mult9(X[i*4+3]);
		d[1] =  Mult9(X[i*4])^Mult14(X[i*4+1])^Mult11(X[i*4+2])^Mult13(X[i*4+3]);
		d[2] = Mult13(X[i*4])^ Mult9(X[i*4+1])^Mult14(X[i*4+2])^Mult11(X[i*4+3]);
		d[3] = Mult11(X[i*4])^Mult13(X[i*4+1])^ Mult9(X[i*4+2])^Mult14(X[i*4+3]);
		
		X[i*4]   = d[0];
		X[i*4+1] = d[1];
		X[i*4+2] = d[2];
		X[i*4+3] = d[3];
	}
} 

void AES_128(byte X[16],byte K[16],byte Y[16]){
	int i;
	byte RK[11][16];
	

	
	memcpy(Y,X,16);
	
	KeyExpansion(K,RK);	
	
	AddRoundKey(Y,RK,0);
	
	for(i = 0; i < 9;i++){
		SubBytes(Y);	
		
		ShiftRows(Y);	
		
		MixColumns(Y);	
		
		AddRoundKey(Y,RK,i+1);
	}
	
	SubBytes(Y);
	
	ShiftRows(Y);
	
	AddRoundKey(Y,RK,i+1);
	
	memcpy(K,RK[0],16);
}

void AES_128_1(byte X[16],byte K[16],byte Y[16]){
	int i;
	byte RK[11][16];
	
	
	
	memcpy(Y,X,16);

	KeyExpansion(K,RK);	
	
	AddRoundKey(Y,RK,10);
	
	for(i = 9; i > 0;i--){
		ShiftRows_1(Y);	
		
		SubBytes_1(Y);	
		
		AddRoundKey(Y,RK,i);
		
		MixColumns_1(Y);	
	}
	
	ShiftRows_1(Y);
	
	SubBytes_1(Y);
	
	AddRoundKey(Y,RK,0);
	
	
	memcpy(K,RK[0],16);
}

void BicliqueKeyExpansion(byte K[16],byte RK[3][16]){
	int i;
	
	
	for(i = 8; i < 10;i++){
		memcpy(RK[i-8],K,16);
		KeySchedule(K,i);
	}
	memcpy(RK[2],K,16);
	memcpy(K,RK[0],16);

}

void BicliqueKeyExpansion_1(byte K[16],byte RK[9][16]){
	int i;
	

	for(i = 7; i >=0;i--){
		memcpy(RK[i+1],K,16);
		KeySchedule_1(K,i);
	}
	memcpy(RK[0],K,16);
	memcpy(K,RK[8],16);
}

void g(byte X[16],byte K[16],byte Y[16]){
	int i;
	byte RK[9][16];
	

	
	memcpy(Y,X,16);
	
	BicliqueKeyExpansion_1(K,RK);	
	
	AddRoundKey(Y,RK,0);
	for(i = 0; i < 7;i++){
		SubBytes(Y);	
		ShiftRows(Y);	
		MixColumns(Y);	
		
		AddRoundKey(Y,RK,i+1);
	}
}

void f(byte X[16],byte K[16],byte Y[16]){
	int i;
	byte RK[3][16];
	

	
	memcpy(Y,X,16);
	
	BicliqueKeyExpansion(K,RK);	
	
	for(i = 0; i < 2;i++){
		SubBytes(Y);	
		
		ShiftRows(Y);	
		
		MixColumns(Y);	
		
		AddRoundKey(Y,RK,i);
	}
	
	SubBytes(Y);
	
	ShiftRows(Y);
	
	AddRoundKey(Y,RK,i);
	
	memcpy(K,RK[0],16);
}

void f_1(byte X[16],byte K[16],byte Y[16]){
	int i;
	byte RK[3][16];
	

	
	memcpy(Y,X,16);
	
	
	BicliqueKeyExpansion(K,RK);	
	
	AddRoundKey(Y,RK,2);
	
	ShiftRows_1(Y);
	
	SubBytes_1(Y);
	
	for(i = 1; i >=0;i--){
		
		AddRoundKey(Y,RK,i);
		
		MixColumns_1(Y);
		
		ShiftRows_1(Y);	
		
		SubBytes_1(Y);	
	}
	memcpy(K,RK[0],16);
} 

byte r(byte X[16],byte precomp[7],byte RK[3][16]){
	int i,j=0;
	byte Y[16];

	memcpy(Y,X,16);

	AddRoundKey(Y,RK,0);

	SubBytes(Y);	
	
	ShiftRows(Y);	
	
	///guardar este estado	
	precomp[0] = Y[1];
	precomp[1] = Y[2];
	precomp[2] = Y[3];
	precomp[3] = Y[6];
	precomp[4] = Y[9];
	precomp[5] = Y[11];
	precomp[6] = Y[15];
	
	MixColumns(Y);	
	
	AddRoundKey(Y,RK,1);

	SubBytes(Y);	
	
	ShiftRows(Y);	
	
	MixColumns(Y);	

	AddRoundKey(Y,RK,2);
	
	return Y[12];
}

byte t_1(byte X[16],byte precomp[2][16],byte RK[9][16]){
	int i,j=0;
	byte a,Y[16];
	
	memcpy(Y,X,16);

	AddRoundKey(Y,RK,7);
	
	for(i = 6; i >=3;i--){
		
		MixColumns_1(Y);
	
		if(i == 6) memcpy(precomp[0],Y,16);
		
		
		ShiftRows_1(Y);	
	
		if(i == 6) memcpy(precomp[1],Y,16);
				
		
		
		SubBytes_1(Y);

		if(i == 6) memcpy(precomp[2],Y,16);
	
		AddRoundKey(Y,RK,i);
		if(i == 6)	memcpy(precomp[3],Y,16);
				
	}
	
	MixColumns_1(Y);
		
	ShiftRows_1(Y);	
	
	SubBytes_1(Y);	
	
	return Y[12];
}

