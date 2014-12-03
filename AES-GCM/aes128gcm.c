
#include "aes128gcm.h"

/*************************************************************
 ********* DECLARATION OF RECOMMENDATION FUNCTIONS ***********
 *************************************************************/

void inc32(unsigned char* incJ, unsigned char* J, const unsigned long len);
void dot(unsigned char* out, unsigned char* A, unsigned char* B);
void ghash(unsigned char* S, unsigned char* uncS, unsigned char* H, const unsigned long len);
void gctr(unsigned char* ret, const unsigned char *k, const unsigned char *ICB, const unsigned char *X, const unsigned long len);

/*************************************************************
 ***** AUXILIARY FUNCTIONS FOR INTERNAL REPRESENTATIONS ******
 *************************************************************/

/*
 Converts a 128 bit char array to a 128*sizeof(unsigned int) int array.
 Each integer in the array represents one bit of the original char array.
 Used in dot.
 */
unsigned int* byteArray(unsigned char* chars) {
    unsigned int *res = (unsigned int*) calloc(sizeof(unsigned int), 128);
    int x = 0;
    int i;
    for(i = 0; i < 16; i++) {
        int k;
        char m = 0x01;
        for (k = 7; k >= 0; k--) {
            char r = ((m << k) & chars[i])>>k;
            res[x++] = r;
        }
    }
    return res;
}

/*
 Converts a 128*sizeof(unsigned int) int array to a 128 bit char array.
 Each integer in the array represents one bit of the new char array.
 Used in dot.
 */
unsigned char* backToChar(unsigned int* byteArray) {
    unsigned char* res = (unsigned char*) calloc(sizeof(unsigned char), 16);
    int pow[] = {1, 2, 4, 8, 16, 32, 64, 128};
    int i;
    for (i = 0; i < 16; ++i) {
        int k;
        for (k = 0; k < 8; ++k) {
            res[i] += byteArray[i*8+k]*pow[7-k];
        }
    }
    return res;
}

/*************************************************************
 **** AUXILIARY FUNCTIONS FOR OPERATIONS WITH INT ARRAYS *****
 *************************************************************/

/*
 Performs a shift right in the int array.
 Used in dot.
 */
void shiftRight(unsigned int* V) {
    int i;
    for (i = 127; i > 0; i--) {
        V[i] = V[i-1];
    }
    V[0] = 0;
}

/*
 Performs a xor on the int arrays and returns the result.
 Used in dot.
 */
unsigned int* xor(unsigned int* a, unsigned int* b) {
    unsigned int* ret = (unsigned int*) calloc(sizeof(unsigned int), 128);
    int i;
    for (i = 0; i < 128; i++) {
        ret[i] = a[i] ^ b[i];
    }
    return ret;
}

/*************************************************************
 * AUXILIARY FUNCTIONS FOR CONVERSIONS BETWEEN INT AND CHAR* *
 *************************************************************/

/*
 Converts a 32 bit char array to an integer.
 Used in inc32.
 */
int toInt(unsigned char* lsb) {
    int pow[] = {1,256,65536,16777216};
    int val = 0;
    int i;
    for (i = 3; i >= 0; i--) {
        val += lsb[3-i] * pow[i];
    }
    return val;
}

/*
 Converts an integer to a 32 bit char array.
 Used in inc32 and to obtain the bit string representation of the lenght of A anc C.
 Converts little endian int representation in x86 to big endian bit string.
 */
void fromInt(unsigned char* toCopy, int val) {
    unsigned char* temp = (unsigned char*) calloc(sizeof(unsigned char), 4);
    memcpy(temp, &val, 4*sizeof(unsigned char));
    int i;
    for (i = 0; i < 4; i++) {
        memcpy(&toCopy[i], &temp[3-i], 1);
    }
    free(temp);
}

/*************************************************************
 ****** AUXILIARY FUNCTIONS FOR PERFORMING COMPUTATIONS ******
 *************************************************************/

/*
 Implements the computation defined in line 3 of 6.4.
 Used in ghash.
 */
void ghash_int(unsigned char* Y1, unsigned char* Y0, unsigned char* S, unsigned char* H) {
    // XOR part of Yi = (Yi-1 XOR Xi) . H
    int i;
    for (i = 0; i < 16; i++) {
        Y1[i] = Y0[i] ^ S[i];
    }
    
    // Multiplication part of Yi = (Yi-1 XOR Xi) . H
    unsigned char* out = (unsigned char*) calloc(sizeof(unsigned char), 16);
    dot(out, Y1, H);
    
    // Return result of computation
    memcpy(Y1, out, 16*sizeof(unsigned char));
    free(out);
}

/*************************************************************
 ************* RECOMMENDATION DEFINED FUNCTIONS **************
 *************************************************************/

/*
 Implements the incrementing function defined in 6.2.
 Used in aes128gcm.
 */
void inc32(unsigned char* incJ, unsigned char* J, const unsigned long len) {
    // Obtains the 32 bit LSB of J
    unsigned char* lsb = (unsigned char*) calloc(sizeof(unsigned char), 4);
    memcpy(lsb, &J[len*16-4], 4);
    
    // Calculates increment mod 2^32
    int lsbInt = toInt(lsb);
    lsbInt++;
    lsbInt = lsbInt % (long)pow(2, 32);
    
    // Concatenates the 96 bit MSB with the increment
    memcpy(incJ, J, len*16-4);
    unsigned char* intt = &incJ[len*16-4];
    fromInt(intt, lsbInt);
    free(lsb);
}

/*
 Implements the multiplication operation on blocks defined in 6.3.
 A = X, B = Y, C = R
 Used in aes128gcm.
 */
void dot(unsigned char* out, unsigned char* A, unsigned char* B) {
    unsigned char* C = (unsigned char*) calloc(sizeof(unsigned char), 16);
    C[0] = 0xE1;
    
    // Convert the char arrays to internal int array representation
    unsigned int* X = byteArray(A);
    unsigned int* Y = byteArray(B);
    unsigned int* R = byteArray(C);
    
    free(C);
    
    // Z0 = 0^128; Y = Y (to keep the same names) (line 2)
    unsigned int* Z = (unsigned int*) calloc(sizeof(unsigned int), 128);
    unsigned int* V = Y;
    
    // For i = 0 ... 127 (line 3)
    int i;
    for (i = 0; i < 128; i++) {
        // Computation of Z
        if(X[i] == 1) {
            unsigned int* t = xor(Z, V);
            memcpy(Z, t, 128*sizeof(unsigned int));
            free(t);
        }
        
        // Computation of V
        if(V[127] == 0) {
            shiftRight(V);
        } else {
            shiftRight(V);
            unsigned int* t = xor(V, R);
            memcpy(V, t, 128*sizeof(unsigned int));
            free(t);
        }
    }
    
    free(X);
    free(Y);
    free(R);
    
    // Convert back to char array
    unsigned char* ret = backToChar(Z);
    
    free(Z);
    
    // Return Z128 (line 4)
    memcpy(out, ret, 16);
    
    free(ret);
}

/*
 Implements the GHASH function defined in 6.4.
 Used in aes128gcm.
 */
void ghash(unsigned char* S, unsigned char* uncS, unsigned char* H, const unsigned long len) {
    // Y0 = 0^128 (assured by the allocation of memory using the calloc function)
    unsigned char* Y = (unsigned char*) calloc(sizeof(unsigned char), (len+1)*16);
    
    // For i = 1 ... m (line 3)
    int i;
    for (i = 0; i < len; i++) {
        // Yi = (Yi-1 XOR Xi) . H (line 3)
        ghash_int(&Y[(i+1)*16], &Y[i*16], &uncS[i*16], H);
    }
    
    // Return Ym (line 4)
    memcpy(S, &Y[len*16], 16);
    free(Y);
}

/*
 Implements the GCTR function defined in 6.5.
 Used in aes128gcm.
 */
void gctr(unsigned char* ret, const unsigned char *k, const unsigned char *ICB, const unsigned char *X, const unsigned long len) {
    if(len != 0) {
        // Allocate space for the CB
        unsigned char* CB = (unsigned char*) calloc(sizeof(unsigned char), 16*len);
        
        // CB1 = ICB
        memcpy(CB, ICB, sizeof(unsigned char)*16);
        
        // For i = 2 ... n (line 5)
        int i;
        for (i = 16; i < len*16; i += 16) {
            inc32(&CB[i], &CB[i-16], 1);
        }
        
        // For i = 1 ... n -1 (line 6)
        for (i = 0; i < len*16; i += 16) {
            unsigned char* temp = (unsigned char*) calloc(sizeof(unsigned char), 16);
            aes128e(temp, &CB[i], k);
            int k;
            for (k = 0; k < 16; k++) {
                temp[k] = temp[k] ^ X[i+k];
            }
            memcpy(&ret[i], temp, 16);
            free(temp);
        }
        free(CB);
    }
}

/*************************************************************
 *********** IMPLEMENTATION OF THE RECOMMENDATION ************
 *************************************************************/

/* 
 Under the 16-byte (128-bit) key "k",
 and the 12-byte (96-bit) initial value "IV",
 encrypt the plaintext "plaintext" and store it at "ciphertext".
 The length of the plaintext is a multiple of 16-byte (128-bit) given by len_p (e.g., len_p = 2 for a 32-byte plaintext).
 The length of the ciphertext "ciphertext" is len_p*16 bytes.
 The authentication tag is obtained by the 16-byte tag "tag".
 For the authentication an additional data "add_data" can be added.
 The number of blocks for this additional data is "len_ad" (e.g., len_ad = 1 for a 16-byte additional data).
 */
void aes128gcm(unsigned char *ciphertext, unsigned char *tag, const unsigned char *k, const unsigned char *IV, const unsigned char *plaintext, const unsigned long len_p, const unsigned char* add_data, const unsigned long len_ad) {
    
    // Calculate H
    unsigned char *temp = (unsigned char*) calloc(sizeof(unsigned char), 16);
    unsigned char *H = (unsigned char*) calloc(sizeof(unsigned char), 16);
    aes128e(H, temp, k);
    
    // Calculate J0
    unsigned char *J = temp;
    memcpy(J, IV, sizeof(unsigned char)*12);
    J[15] = 1;
    
    // Calculate C,
    unsigned char *incJ = (unsigned char*) calloc(sizeof(unsigned char), 16);
    inc32(incJ, J, 1);
    gctr(ciphertext, k, incJ, plaintext, len_p);
    
    // u = 0 and v = 0 (since the size of the cyphertext and of the additional data is a multiple of 128, they are always 0)
    
    // Calculate S
    unsigned char *uncS = (unsigned char*) calloc(sizeof(unsigned char), (1+len_ad+len_p)*16);
    memcpy(uncS, add_data, 16*len_ad);
    memcpy(&uncS[16*len_ad], ciphertext, 16*len_p);
    
    unsigned int lenA = (unsigned int)len_ad*128;
    unsigned int lenC = (unsigned int)len_p*128;
    
    fromInt(&uncS[(len_ad+len_p)*16+4], lenA);
    fromInt(&uncS[(len_ad+len_p)*16+12], lenC);
    
    unsigned char *S = (unsigned char*) calloc(sizeof(unsigned char), 16);
    ghash(S, uncS, H, 1+len_ad+len_p);
    
    // Calculate T
    gctr(tag, k, J, S, 1);
    
    free(temp);
    free(H);
    free(incJ);
    free(uncS);
    free(S);
}