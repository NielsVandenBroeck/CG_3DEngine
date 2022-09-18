//
// Created by Niels on 26/02/2021.
//

#ifndef ENGINE_L2D_H
#define ENGINE_L2D_H

#include"Parsing/l_parser.h"
#include "Line2D.h"
#include "Libraries/easy_image.h"
#include<stack>
#include<cmath>


class L2D {
public:
    L2D(const LParser::LSystem2D &l_system);
    vector<Line2D*> startIteration();
    void setLineColor(const vector<double> &linecolor);
    vector<double> getMinMax();
private:
    vector<Line2D*> iterate(const string iteration,const int nr_iteration);
    double degree_to_rad(const double degree);
    Line2D* draw();
    void skip();

    double angle;
    set<char> alphabet;
    string initiator;
    map<char,string> rules;
    map<char,int> draws;
    unsigned int nr_iterations;
    vector<double> lineColor;

    double currentAngle = 0;
    pair<double,double> CurrentPosition;

    stack<double> angleStack;
    stack<pair<double,double>> positionStack;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    bool first = true;

};


#endif //ENGINE_L2D_H
