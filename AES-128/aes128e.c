#include <stdint.h>
#include "aes128e.h"
#include <stdlib.h>
#include <string.h>

/* Multiplication by two in GF(2^8). Multiplication by three is xtime(a) ^ a */
#define xtime(a) ( ((a) & 0x80) ? (((a) << 1) ^ 0x1b) : ((a) << 1) )

// Macro that provides multiplication by three as defined above
#define txtime(a) (xtime(a)^a)

// Number of rounds defined by the standard for AES-128
#define NROUNDS 10

/* The S-box table */
static const unsigned char sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

/* The round constant table (needed in KeyExpansion) */
static const unsigned char rcon[10] = {
    0x01, 0x02, 0x04, 0x08, 0x10,
    0x20, 0x40, 0x80, 0x1b, 0x36 };

// This function implements the bitwise XOR operation between two char arrays of lenght 4
void bitwise_xor(unsigned char* temp, unsigned char* to_xor) {
    int i;
    for (i = 0; i < 4; i++) {
        temp[i] = temp[i]^to_xor[i];
    }
}

// This function implements the bitwise XOR operation between two char arrays of arbitrary lenght
// where one of them is constant (implemented to eliminate warnings related to losing classifiers)
void bitwise_xor_const(unsigned char* temp, const unsigned char* to_xor, int size) {
    int i;
    for (i = 0; i < size; i++) {
        temp[i] = temp[i]^to_xor[i];
    }
}

// Implements the SubBytes() transformation from section 5.1.1 of the AES standard
// using the provided sbox array
unsigned char sub_byte(unsigned char byte) {
    return sbox[byte];
}

// Implements the AddRoundKey() transformation from section 5.1.4 of the AES standard
//
// Due to the internal representation of the array of round keys, the column must be
// transferred to a temporary array in order to XOR with the round key.
void add_round_key_col(unsigned char *state, const unsigned char *key, int col) {
    unsigned char temp[4];
    temp[0] = state[col];
    temp[1] = state[col+4];
    temp[2] = state[col+8];
    temp[3] = state[col+12];
    bitwise_xor_const(temp, key, 4);
    state[col] = temp[0];
    state[col+4] = temp[1];
    state[col+8] = temp[2];
    state[col+12] = temp[3];
}

// Applies the AddRoundKey() transformation to all the columns in the state
void add_round_key(unsigned char *state, const unsigned char *key_schedule) {
    int i;
    for (i=0; i<4; i++) {
        add_round_key_col(state, &key_schedule[i*4], i);
    }
}

// Applies the SubBytes() transformation to all the bytes in the state
void sub_bytes(unsigned char *state) {
    int i;
    for (i = 0; i < 16; i++) {
        state[i] = sub_byte(state[i]);
    }
}

// Implements the ShiftRows() transformation from section 5.1.2 of the AES standard
void shift_rows(unsigned char *state) {
    
    // To keep the previous state while we are shifting the current one, allowing to
    // correctly shift it, the state is copied at the beginning
    unsigned char copy_state[16];
    memcpy(copy_state, state, 16);
    
    // Row 1 shift
    state[4] = copy_state[5];
    state[5] = copy_state[6];
    state[6] = copy_state[7];
    state[7] = copy_state[4];
    
    // Row 2 shift
    state[8] = copy_state[10];
    state[9] = copy_state[11];
    state[10] = copy_state[8];
    state[11] = copy_state[9];
    
    // Row 3 shift
    state[12] = copy_state[15];
    state[13] = copy_state[12];
    state[14] = copy_state[13];
    state[15] = copy_state[14];
}

// Implements the MixColumns() transformation from section 5.1.3 of the AES standard
void mix_column(unsigned char *state, int column) {
    
    // To keep the starting values of the column to perform the calculations, they
    // are copied to temporary variables
    char s0 = state[column];
    char s1 = state[4+column];
    char s2 = state[8+column];
    char s3 = state[12+column];
    
    state[column] = xtime(s0)^txtime(s1)^s2^s3;
    state[4+column] = s0^xtime(s1)^txtime(s2)^s3;
    state[8+column] = s0^s1^xtime(s2)^txtime(s3);
    state[12+column] = txtime(s0)^s1^s2^xtime(s3);
}

// Applies the MixColumns() transformation to all the columns in the state
void mix_columns(unsigned char *state) {
    int i;
    for (i = 0; i < 4; i++) {
        mix_column(state, i);
    }
}

// Implements the SubWord() transformation from section 5.2 of the AES standard
// Applies the SubBytes() transformation to all the bytes in a word
void sub_word(unsigned char* word) {
    int i;
    for (i = 0; i < 4; i++) {
        word[i] = sub_byte(word[i]);
    }
}

// Implements the RotWord() transformation from section 5.2 of the AES standard
void rot_word(unsigned char* word) {
    char temp = word[0];
    word[0] = word[1];
    word[1] = word[2];
    word[2] = word[3];
    word[3] = temp;
}

// Implements the key expansion described in section 5.2 of the AES standard
void key_expansion(const unsigned char *key, unsigned char *key_schedule) {
    
    unsigned char temp[4];
    
    // Instead of performing the first while cycle, the key is copied to the
    // beginning of the key_schedule, which produces the same result
    memcpy(key_schedule, key, 16);
    
    int i = 4;
    while (i < 44) {
        memcpy(temp, &key_schedule[(i-1)*4], 4);
        if (i % 4 == 0) {
            rot_word(temp);
            sub_word(temp);
            bitwise_xor_const(temp, &rcon[i/4-1],1);
        }
        bitwise_xor(temp, &key_schedule[(i-4)*4]);
        memcpy(&key_schedule[i*4], temp, 4);
        i++;
    }
}

// Copies the input to the state and adapts it to the state representation of AES
void copy_to_state(const unsigned char *p, unsigned char *state) {
    int i;
    for (i = 0; i < 4; i++) {
        state[0+i] = p[4*i+0];
        state[4+i] = p[4*i+1];
        state[8+i] = p[4*i+2];
        state[12+i] = p[4*i+3];
    }
}

// Copies the state to the output and adapts it to the expected array representation
void copy_from_state(unsigned char *state, unsigned char *c) {
    int i;
    for (i = 0; i < 4; i++) {
        c[0+i*4] = state[0+i];
        c[1+i*4] = state[4+i];
        c[2+i*4] = state[8+i];
        c[3+i*4] = state[12+i];
    }
}

// Implements the cipher in section 5.1 of the AES standard
/* Under the 16-byte key at k, encrypt the 16-byte plaintext at p and store it at c. */
void aes128e(unsigned char *c, const unsigned char *p, const unsigned char *k) {
    
    // Key expansion
    unsigned char* key_schedule = (unsigned char*) calloc(176,1);
    key_expansion(k, key_schedule);
    
    // Initialization of the AES state
    unsigned char state[16];
    copy_to_state(p,state);
    
    // "Round 0" of the AES-128
    add_round_key(state, key_schedule);
    
    // Rounds 1 to 9 of the AES-128
    int rounds;
    for (rounds = 1; rounds < NROUNDS; rounds++) {
        sub_bytes(state);
        shift_rows(state);
        mix_columns(state);
        add_round_key(state, &key_schedule[rounds*16]);
    }
    
    // Round 10 of the AES-128
    sub_bytes(state);
    shift_rows(state);
    add_round_key(state, &key_schedule[160]);
    
    free(key_schedule);
    
    // Copies the state to the output variable in the correct format
    copy_from_state(state, c);
}

