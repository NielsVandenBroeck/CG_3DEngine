//
// Created by Niels on 16/03/2021.
//

#define _USE_MATH_DEFINES
#include "L3D.h"

L3D::L3D(const LParser::LSystem3D &l_system) {
    alphabet = l_system.get_alphabet();
    for(auto i = alphabet.begin(); i != alphabet.end(); i++){
        draws[*i] = l_system.draw(*i);
        rules[*i] = l_system.get_replacement(*i);
    }
    initiator = l_system.get_initiator();
    angleInRad = l_system.get_angle() * (M_PI / 180);
    nr_iterations = l_system.get_nr_iterations();
}

img::Figure* L3D::startIteration() {
    return iterate(initiator, nr_iterations);
}

img::Figure* L3D::iterate(const string iterator, const int nr_iteration) {
    if(nr_iteration == nr_iterations){
        f->points.push_back(CurrentPosition);
        f->faces.push_back({0});
    }
    for(int i = 0; i < iterator.size(); i++){
        if(iterator[i] == '+'){
            Vector3D h = currentH;
            Vector3D l = currentL;
            currentH = {h*cos(angleInRad)+l*sin(angleInRad)};
            currentL = {-h*sin(angleInRad)+l*cos(angleInRad)};
        }
        else if(iterator[i] == '-') {
            Vector3D h = currentH;
            Vector3D l = currentL;
            currentH = {h*cos(-angleInRad)+l*sin(-angleInRad)};
            currentL = {-h*sin(-angleInRad)+l*cos(-angleInRad)};
        }
        else if(iterator[i] == '^'){
            Vector3D h = currentH;
            Vector3D u = currentU;
            currentH = {h*cos(angleInRad)+u*sin(angleInRad)};
            currentU = {-h*sin(angleInRad)+u*cos(angleInRad)};
        }
        else if(iterator[i] == '&') {
            Vector3D h = currentH;
            Vector3D u = currentU;
            currentH = {h*cos(-angleInRad)+u*sin(-angleInRad)};
            currentU = {-h*sin(-angleInRad)+u*cos(-angleInRad)};
        }
        else if(iterator[i] == '\\'){
            Vector3D l = currentL;
            Vector3D u = currentU;
            currentL = {l*cos(angleInRad)-u*sin(angleInRad)};
            currentU = {l*sin(angleInRad)+u*cos(angleInRad)};
        }
        else if(iterator[i] == '/') {
            Vector3D l = currentL;
            Vector3D u = currentU;
            currentL = {l*cos(-angleInRad)-u*sin(-angleInRad)};
            currentU = {l*sin(-angleInRad)+u*cos(-angleInRad)};
        }
        else if(iterator[i] == '|') {
            currentH = -currentH;
            currentL= -currentL;
        }
        else if(iterator[i] == '('){
            angleStack.push({currentH,currentL,currentU});
            positionStack.push(CurrentPosition);
        }
        else if(iterator[i] == ')'){
            currentH = angleStack.top()[0];
            currentL = angleStack.top()[1];
            currentU = angleStack.top()[2];
            angleStack.pop();
            CurrentPosition = positionStack.top();
            f->points.push_back(CurrentPosition);
            int size = f->points.size();
            f->faces.push_back({size-1});
            faceindex++;
            positionStack.pop();
        }
        else if(nr_iteration > 0){
            iterate(rules[iterator[i]],nr_iteration-1);
        }
        else{
            set<char>::iterator it = alphabet.begin();
            while(true){
                if(*it == iterator[i]){
                    bool foundnext = false;
                    if(i < iterator.size()-1){
                        set<char>::iterator it2 = alphabet.begin();
                        while(!foundnext && it2 != alphabet.end()) {
                            if (*it2 == iterator[i + 1]) {
                                foundnext = true;
                                break;
                            }
                            it2++;
                        }
                    }
                    if(foundnext){
                        nextPos();
                    }
                    else{
                        f->points.push_back(nextPos());
                        if(draws[*it] == 1){
                            int index = f->points.size()-1;
                            f->faces[faceindex].push_back(index);
                        }
                    }
                    break;
                }
                it++;
            }
        }
    }
    return f;
}

Vector3D L3D::nextPos() {
    CurrentPosition = Vector3D::point(CurrentPosition+currentH);
    return CurrentPosition;
}
