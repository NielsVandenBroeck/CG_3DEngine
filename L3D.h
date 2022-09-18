//
// Created by Niels on 16/03/2021.
//

#ifndef ENGINE_L3D_H
#define ENGINE_L3D_H

#include"Libraries/easy_image.h"
#include"Parsing/l_parser.h"
#include"Libraries/vector3d.h"
#include<vector>
#include<cmath>
#include<stack>


using namespace std;

class L3D {
public:
    L3D(const LParser::LSystem3D &l_system);
    img::Figure* startIteration();
    img::Figure* iterate(const string iteration,const int nr_iteration);

    Vector3D nextPos();

private:
    double angleInRad;

    Vector3D currentH = Vector3D::vector(1,0,0);
    Vector3D currentL = Vector3D::vector(0,1,0);
    Vector3D currentU = Vector3D::vector(0,0,1);

    int faceindex = 0;

    img::Figure* f = new img::Figure();
    set<char> alphabet;
    string initiator;
    map<char,string> rules;
    map<char,int> draws;
    unsigned int nr_iterations;

    Vector3D CurrentPosition = Vector3D::point(0,0,0);

    stack<vector<Vector3D>> angleStack;
    stack<Vector3D> positionStack;
};


#endif //ENGINE_L3D_H
