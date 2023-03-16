#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

uint8_t sbox[256] = {
        0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
        0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
        0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
        0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
        0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
        0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
        0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
        0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
        0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
        0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
        0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
        0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
        0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
        0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
        0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
        0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};


const uint8_t rcon[10] = { 1, 2, 4, 8, 16, 32, 64, 128, 27, 54 };
uint8_t key[16];

uint8_t masked_sbox[256];
uint8_t m1;
uint8_t m2;
uint8_t mask1[4];
uint8_t mask2[4]; // mask2 = mix_column(mask1)


void display_vector(uint8_t v[16]) {
    for (int32_t i = 0; i < 16; i += 1) {
        printf("%.2x", v[i]);
    }
    printf("\n");
}

void display_x(uint8_t x[4][4]) {
    printf("\n");
    for (int32_t i = 0; i < 16; i += 1) {
        printf("%.2x", x[i % 4][i / 4]);
    }
    printf("\n");
}


void shift_rows(uint8_t x[4][4]) {
    uint8_t tmp = x[1][0];
    x[1][0] = x[1][1];
    x[1][1] = x[1][2];
    x[1][2] = x[1][3];
    x[1][3] = tmp;

    tmp = x[2][0];
    x[2][0] = x[2][2];
    x[2][2] = tmp;;
    tmp = x[2][1];
    x[2][1] = x[2][3];
    x[2][3] = tmp;

    tmp = x[3][0];
    x[3][0] = x[3][3];
    x[3][3] = x[3][2];
    x[3][2] = x[3][1];
    x[3][1] = tmp;
}


/* Copied from wikipedia */
void mix_column(uint8_t r[4]) {
    uint8_t a[4];
    uint8_t b[4];
    /* The array 'a' is simply a copy of the input array 'r'
     * The array 'b' is each element of the array 'a' multiplied by 2
     * in Rijndael's Galois field
     * a[n] ^ b[n] is element n multiplied by 3 in Rijndael's Galois field */
    for (int8_t c = 0; c < 4; c++) {
        a[c] = r[c];
        /* h is 0xff if the high bit of r[c] is set, 0 otherwise */
        uint8_t h = (uint8_t) ((signed char) r[c] >> 7); /* arithmetic right shift, thus shifting in either zeros or ones */
        b[c] = r[c] << 1; /* implicitly removes high bit because b[c] is an 8-bit char, so we xor by 0x1b and not 0x11b in the next line */
        b[c] ^= 0x1B & h; /* Rijndael's Galois field */
    }
    r[0] = b[0] ^ a[3] ^ a[2] ^ b[1] ^ a[1]; /* 2 * a0 + a3 + a2 + 3 * a1 */
    r[1] = b[1] ^ a[0] ^ a[3] ^ b[2] ^ a[2]; /* 2 * a1 + a0 + a3 + 3 * a2 */
    r[2] = b[2] ^ a[1] ^ a[0] ^ b[3] ^ a[3]; /* 2 * a2 + a1 + a0 + 3 * a3 */
    r[3] = b[3] ^ a[2] ^ a[1] ^ b[0] ^ a[0]; /* 2 * a3 + a2 + a1 + 3 * a0 */
}


void add_round_key(uint8_t x[4][4], int32_t round, uint8_t round_key[176]) {
    for (int32_t i = 0; i < 4; i += 1) {
        for (int32_t j = 0; j < 4; j += 1) {
            x[j][i] ^= round_key[round * 16 + i * 4 + j];
        }
    }
}


void key_schedule(uint8_t key[16], uint8_t round_key[176]) {
    /* Init round_key by master key */
    for (int32_t i = 0; i < 4; i += 1) {
        round_key[i * 4]     = key[i * 4];
        round_key[i * 4 + 1] = key[i * 4 + 1];
        round_key[i * 4 + 2] = key[i * 4 + 2];
        round_key[i * 4 + 3] = key[i * 4 + 3];
    }
    for (int32_t i = 4; i < 44; i += 1) {
        uint8_t temp[4];
        for (int32_t j = 0; j < 4; j += 1) {
            temp[j] = round_key[(i - 1) * 4 + j];
        }
        if (i % 4 == 0) {
            /* RotWord */
            uint8_t k = temp[0];
            temp[0] = temp[1];
            temp[1] = temp[2];
            temp[2] = temp[3];
            temp[3] = k;

            /* S-box */
            temp[0] = sbox[temp[0]];
            temp[1] = sbox[temp[1]];
            temp[2] = sbox[temp[2]];
            temp[3] = sbox[temp[3]];

            /* RCON */
            temp[0] = temp[0] ^ rcon[(i / 4) - 1];
        }
        round_key[i * 4 + 0] = round_key[(i - 4) * 4 + 0] ^ temp[0];
        round_key[i * 4 + 1] = round_key[(i - 4) * 4 + 1] ^ temp[1];
        round_key[i * 4 + 2] = round_key[(i - 4) * 4 + 2] ^ temp[2];
        round_key[i * 4 + 3] = round_key[(i - 4) * 4 + 3] ^ temp[3];
    }

}

void mask_key(uint8_t round_key[176]){
    for(int32_t round = 0; round < 11; ++round){
        for(int32_t i = 0; i < 4; ++i){
            for(int32_t j = 0; j < 4; ++j){
                if (round == 10){
                    round_key[round * 16 + i * 4 + j] ^= m2;
                } else {
                    round_key[round * 16 + i * 4 + j] ^= mask2[j] ^ m1;
                }
            }
        }
    }
}

void sub_byte(uint8_t x[4][4], uint8_t sbox[256]) {
    for (int32_t i = 0; i < 4; i += 1) {
        for (int32_t j = 0; j < 4; j += 1) {
            x[i][j] = sbox[x[i][j]];
        }
    }
}

void sboxmask_generate(void ){
    for(int32_t i = 0; i < 256; ++i){
        masked_sbox[i ^ m1] = sbox[i] ^ m2; // s'[t ^ m] = s[t] ^ m' , d'après l'article
    }
}

void changemask(uint8_t mask[4][4]){
    for (int32_t i = 0; i < 4; ++i) {
        for (int32_t j = 0; j < 4; ++j) {
            mask[i][j] = mask[i][j] ^ mask1[i] ^ m2;
        }
    }
}

void init(uint8_t x[4][4], uint8_t plain_text[16]) {
    for (int32_t i = 0; i < 4; i += 1) {
        for (int32_t j = 0; j < 4; j += 1) {
            x[j][i] = plain_text[i * 4 + j];
        }
    }
}

void init_mask(void){
    m1 = rand() % 0xFF;
    m2 = rand() % 0xFF;
    for (int32_t i = 0; i < 4; ++i) {
        mask1[i] = rand() % 0xFF;
        mask2[i] = mask1[i];
    }
}

void aes_no_ks(uint8_t x[4][4], uint8_t round_key[176], uint8_t ciphered_text[16]) {
    //1-on initialise le masque
    init_mask();
    //2- mixcolonne sur le masque
    mix_column(mask2);
    //3- on genère la nouvelle sbox
    sboxmask_generate();
    //4- on masque la clé
    mask_key(round_key);
    //5- on masque le message avant d'etre chiffrer
    for(int32_t i = 0; i < 4; i++){
        for(int32_t j = 0; j < 4; j++){
            x[i][j] ^= mask2[i];
        }
    }

    // Rounds
    for (int32_t round = 0; round < 9; round++) {
        add_round_key(x, round, round_key);
        if ( round == 0){
            printf("ici : \t");
            display_x(x);
        }
        sub_byte(x, masked_sbox);// avec la nouvelle sbox
        shift_rows(x);
        changemask(x);
        {
            uint8_t col[4];
            for (int8_t i = 0; i < 4; i += 1) {
                col[0] = x[0][i];
                col[1] = x[1][i];
                col[2] = x[2][i];
                col[3] = x[3][i];
                mix_column(col);
                x[0][i] = col[0];
                x[1][i] = col[1];
                x[2][i] = col[2];
                x[3][i] = col[3];
            }
        }
    }

    // Final Round (no Mix Columns)
    add_round_key(x, 9, round_key);
    sub_byte(x, masked_sbox); //avec la nouvelle sbox
    shift_rows(x);
    add_round_key(x, 10, round_key);

    // Writing output ciphered text
    for (int32_t i = 0; i < 16; i += 1) {
        ciphered_text[i] = x[i % 4][i / 4];
    }
}


uint8_t get_key(uint8_t *k, uint8_t len) {
    for (int32_t i = 0; i < len; i += 1) {
        key[i] = k[i];
    }
    return 0;
}


uint8_t reset(uint8_t *x, uint8_t len) {
    return 0;
}


uint8_t get_pt(uint8_t* pt, uint8_t len) {
    uint8_t round_key[176];
    uint8_t ciphered_text[16];
    uint8_t x[4][4]; // state
    init(x, pt);
    key_schedule(key, round_key);
    aes_no_ks(x, round_key, ciphered_text);

    display_vector(ciphered_text);

    return 0;
}

int main(void) {
    srand(7);
    uint8_t pt[16] = { 0x77, 0x77, 0x77, 0x77, 0x77, 0x77, 0x77, 0x77, 0x77,
                       0x77, 0x77, 0x77, 0x77, 0x77, 0x77, 0x77};
    get_key(pt, 16); // key
    printf("cipher text : ");
    get_pt(pt, 16);
    printf("cipher text : 0afe465c4e9f30901de565641e77cf7c");
    return 0;
}