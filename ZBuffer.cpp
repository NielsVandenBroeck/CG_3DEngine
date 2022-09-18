//
// Created by Niels on 20/03/2021.
//

#include "ZBuffer.h"
#include <iostream>

ZBuffer::ZBuffer(const int width, const int height) {
    double inf = numeric_limits<double>::infinity();
    vector<vector<double>> b(width, vector<double>(height,inf));
    zBuffer = b;
}