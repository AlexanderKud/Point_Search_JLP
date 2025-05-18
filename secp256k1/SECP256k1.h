#ifndef SECP256K1H
#define SECP256K1H

#include "Point.h"
#include <string>

class Secp256K1 {

public:

  Secp256K1();
  ~Secp256K1();
  void Init();
  Point ScalarMultiplication(Int *privKey);
  Point PointMultiplication(Point &P, Int *scalar);
  Point PointDivision(Point &P, Int *scalar);

  bool  EC(Point &p);

  std::string GetPublicKeyHex(Point &p);
  Point ParsePublicKeyHex(std::string pub);

  
  Point Add(Point &p1, Point &p2);
  Point Add2(Point &p1, Point &p2);
  Point AddPoints(Point &p1, Point &p2);
  Point Double(Point &p);
  Point DoublePoint(Point &p);
  Point Subtract(Point &p1, Point &p2);
  Point SubtractPoints(Point &p1, Point &p2);

  Point G;                 // Generator
  Int   order;             // Curve order
  Int P;                   // Field Order

private:

  uint8_t GetByte(std::string &str,int idx);

  Int GetY(Int x, bool isEven);
  Point GTable[256*32];       // Generator table

};

#endif // SECP256K1H
