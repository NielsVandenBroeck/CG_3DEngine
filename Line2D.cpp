//
// Created by Niels on 19/02/2021.
//
#include <cmath>
#include "Line2D.h"

void Line2D::roundPoints() {
    Point1.first = round(Point1.first);
    Point1.second = round(Point1.second);
    Point2.first = round(Point2.first);
    Point2.second = round(Point2.second);
}

void Line2D::addPosition(const double dx, const double dy) {
    Point1.first = Point1.first + dx;
    Point1.second = Point1.second + dy;
    Point2.first = Point2.first + dx;
    Point2.second = Point2.second + dy;
}

void Line2D::multiplyPoints(double d){
    Point1.first = Point1.first*d;
    Point1.second = Point1.second*d;
    Point2.first = Point2.first*d;
    Point2.second = Point2.second*d;
}

const pair<double, double> &Line2D::getPoint1() const {
    return Point1;
}

const pair<double, double> &Line2D::getPoint2() const {
    return Point2;
}

const img::Color &Line2D::getColor() const {
    return color;
}

void Line2D::setPoint1(const pair<double, double> &point1) {
    Point1 = point1;
}

void Line2D::setPoint2(const pair<double, double> &point2) {
    Point2 = point2;
}

void Line2D::setColor(const img::Color &color) {
    Line2D::color = color;
}

double Line2D::getZ1() const {
    return z1;
}

void Line2D::setZ1(double z1) {
    Line2D::z1 = z1;
}

double Line2D::getZ2() const {
    return z2;
}

void Line2D::setZ2(double z2) {
    Line2D::z2 = z2;
}
