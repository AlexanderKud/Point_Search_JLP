#ifndef SECP256K1H
#define SECP256K1H

#include "Point.h"
#include <string>
#include <vector>

// Address type
#define P2PKH  0
#define P2SH   1
#define BECH32 2

class Secp256K1 {

public:

  Secp256K1();
  ~Secp256K1();
  void Init();
  Point ComputePublicKey(Int *privKey);             // Fastest Method JLP Original

  Point ScalarMultiplication(Int *scalar);          // Simple Scalar Multiplication
  Point ScalarMultiplication2(Int *scalar);

  Point FastScalarMultiplication(Int *scalar);      // Fast Scalar Multiplication
  Point FastScalarMultiplication2(Int *scalar);

  Point PointMultiplication(Point &P, Int *scalar); // Point Multiplication
  Point PointMultiplication2(Point &P, Int *scalar);

  Point PointDivision(Point &P, Int *scalar);       // Point Division
  Point PointDivision2(Point &P, Int *scalar);

  Point NextKey(Point &key);
  void Check();
  bool  EC(Point &p);

  void GetHash160(int type,bool compressed,
    Point &k0, Point &k1, Point &k2, Point &k3,
    uint8_t *h0, uint8_t *h1, uint8_t *h2, uint8_t *h3);

  void GetHash160(int type,bool compressed, Point &pubKey, unsigned char *hash);

  std::string GetAddress(int type, bool compressed, Point &pubKey);
  std::string GetAddress(int type, bool compressed, unsigned char *hash160);
  std::vector<std::string> GetAddress(int type, bool compressed, unsigned char *h1, unsigned char *h2, unsigned char *h3, unsigned char *h4);
  std::string GetPrivAddress(bool compressed, Int &privKey );
  std::string GetPublicKeyHex(bool compressed, Point &p);
  Point ParsePublicKeyHex(std::string str, bool &isCompressed);
  Point ParsePublicKeyHex2(std::string pub);

  bool CheckPudAddress(std::string address);

  static Int DecodePrivateKey(char *key,bool *compressed);
  
  Point Add(Point &p1, Point &p2);
  Point Add2(Point &p1, Point &p2);
  Point AddDirect(Point &p1, Point &p2);
  Point AddDirect2(Point &p1, Point &p2);
  Point Double(Point &p);
  Point DoubleDirect(Point &p);
  Point Negation(Point &p);
  Point Subtract(Point &p1, Point &p2);

  Point G;                 // Generator
  Int   order;             // Curve order
  Int P;                   // Field Order
  Point Infinity;          // Point At Infinity

private:

  uint8_t GetByte(std::string &str,int idx);

  Int GetY(Int x, bool isEven);
  Point GTable[256*32];       // Generator table
  Point addend_points[256];

};

#endif // SECP256K1H
