//
// Created by Niels on 26/02/2021.
//

#include "L2D.h"
L2D::L2D(const LParser::LSystem2D &l_system) {
    currentAngle = l_system.get_starting_angle();
    angle = l_system.get_angle();
    alphabet = l_system.get_alphabet();
    set<char>::iterator it = alphabet.begin();
    for(auto i = alphabet.begin(); i != alphabet.end(); i++){
        draws[*i] = l_system.draw(*i);
        rules[*i] = l_system.get_replacement(*i);
    }
    initiator = l_system.get_initiator();
    nr_iterations = l_system.get_nr_iterations();
}

vector<Line2D*> L2D::startIteration() {
    return iterate(initiator, nr_iterations);
}


vector<Line2D*> L2D::iterate(const string iterator,const int nr_iteration) {
    vector<Line2D*> lines;
    for(int i = 0; i<iterator.size(); i++){
        if(iterator[i] == '+'){
            currentAngle += angle;
        }
        else if(iterator[i] == '-') {
            currentAngle -= angle;
        }
        else if(iterator[i] == '('){
            angleStack.push(currentAngle);
            positionStack.push(CurrentPosition);
        }
        else if(iterator[i] == ')'){
            currentAngle = angleStack.top();
            angleStack.pop();
            CurrentPosition = positionStack.top();
            positionStack.pop();
        }
        else if(nr_iteration > 0){
            vector<Line2D*> lines1 = iterate(rules[iterator[i]],nr_iteration-1);
            lines.insert(lines.end(), lines1.begin(), lines1.end());
        }
        else{
            set<char>::iterator it = alphabet.begin();
            while(true){
                if(*it == iterator[i]){
                    if(draws[*it] == 1){
                        lines.emplace_back(draw());
                    }
                    else{
                        skip();
                    }
                    break;
                }
                it++;
            }
        }
    }
    return lines;
}

Line2D* L2D::draw() {

    double cosV = cos(degree_to_rad(currentAngle));
    double sinV = sin(degree_to_rad(currentAngle));

    double x = CurrentPosition.first + cosV;
    double y = CurrentPosition.second + sinV;

    Line2D* l = new Line2D();
    pair<double,double> p1 = make_pair(CurrentPosition.first,CurrentPosition.second);
    pair<double,double> p2 = make_pair(x,y);

    if(first == true){
        if(x<=CurrentPosition.first){
            xmin = x;
            xmax = CurrentPosition.first;
        }
        else{
            xmin = CurrentPosition.first;
            xmax = x;
        }
        if(y<=CurrentPosition.second){
            ymin = y;
            ymax = CurrentPosition.second;
        }
        else{
            ymin = CurrentPosition.second;
            ymax = y;
        }
        first = false;
    }
    else{
        if(CurrentPosition.first<xmin) xmin = CurrentPosition.first;
        else if(CurrentPosition.first>xmax) xmax = CurrentPosition.first;
        if(CurrentPosition.second<ymin) ymin = CurrentPosition.second;
        else if(CurrentPosition.second>ymax) ymax = CurrentPosition.second;

        if(x<xmin) xmin = x;
        else if(x>xmax) xmax = x;
        if(y<ymin) ymin = y;
        else if(y>ymax) ymax = y;
    }

    l->setPoint1(p1);
    l->setPoint2(p2);
    img::Color c(lineColor[0]*255,lineColor[1]*255, lineColor[2]*255);
    l->setColor(c);
    CurrentPosition.first=p2.first;
    CurrentPosition.second=p2.second;
    return l;
}

void L2D::skip() {
    double cosV = cos(degree_to_rad(currentAngle));
    double sinV = sin(degree_to_rad(currentAngle));

    double x = CurrentPosition.first + cosV;
    double y = CurrentPosition.second + sinV;

    pair<double,double> p2 = make_pair(x,y);

    CurrentPosition.first = p2.first;
    CurrentPosition.second = p2.second;
    return;
}

double L2D::degree_to_rad(const double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}

void L2D::setLineColor(const vector<double> &linecolor) {
    lineColor = linecolor;
    return;
}

vector<double> L2D::getMinMax() {
    vector<double> minmax;
    minmax.push_back(xmin);
    minmax.push_back(xmax);
    minmax.push_back(ymin);
    minmax.push_back(ymax);
    return minmax;
}