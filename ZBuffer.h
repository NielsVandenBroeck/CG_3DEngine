//
// Created by Niels on 20/03/2021.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H
#include <vector>
#include <limits>

using namespace std;

class ZBuffer{
public:
    ZBuffer(const int width, const int height);
    vector<vector<double>> zBuffer;
};


#endif //ENGINE_ZBUFFER_H
