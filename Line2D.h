//
// Created by Niels on 19/02/2021.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H

#include<iostream>
#include "Libraries/easy_image.h"
using namespace std;

class Line2D {
public:
    void roundPoints();

    void addPosition(const double dx,const double dy);

    void multiplyPoints(const double d);

    const pair<double, double> &getPoint1() const;

    const pair<double, double> &getPoint2() const;

    const img::Color &getColor() const;

    void setPoint1(const pair<double, double> &point1);

    void setPoint2(const pair<double, double> &point2);

    void setColor(const img::Color &color);

    double getZ1() const;

    void setZ1(double z1);

    double getZ2() const;

    void setZ2(double z2);

private:
    pair<double,double> Point1;
    pair<double,double> Point2;
    img::Color color;
    double z1;
    double z2;
};


#endif //ENGINE_LINE2D_H
