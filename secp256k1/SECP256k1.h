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
  bool  EC(Point &p);

  std::string GetPublicKeyHex(Point &p);
  Point ParsePublicKeyHex(std::string pub);

  
  Point Add(Point &p1, Point &p2);
  Point Add2(Point &p1, Point &p2);
  Point AddDirect(Point &p1, Point &p2);
  Point AddDirect2(Point &p1, Point &p2);
  Point Double(Point &p);
  Point DoubleDirect(Point &p);
  Point Negation(Point &p);
  Point Subtract(Point &p1, Point &p2);
  Point SubtractDirect(Point &p1, Point &p2);

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
