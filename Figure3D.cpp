//
// Created by Niels on 5/03/2021.
//

#include "Figure3D.h"

Figure3D::Figure3D(const vector<Vector3D> &p, const vector<vector<int>> &f, const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular, const double reflection){
    points = p;
    faces = f;
    ambientReflection = ambient;
    diffuseReflection = diffuse;
    specularReflection = specular;
    reflectionCoefficient = reflection;
}

void Figure3D::scale(const double s){
    Matrix scale;
    scale(1,1) = s;
    scale(2,2) = s;
    scale(3,3) = s;
    scale(1, 4) = 0;
    scale(2, 4) = 0;
    scale(3, 4) = 0;
    scale(4, 4) = 1;
    for(int i = 0; i<points.size(); i++){
        points[i] = points[i]*scale;
    }
}

void Figure3D::translate(Vector3D &c){
    Matrix center;
    center(4,1) = c.x;
    center(4,2) = c.y;
    center(4,3) = c.z;
    center(1, 4) = 0;
    center(2, 4) = 0;
    center(3, 4) = 0;
    center(4, 4) = 1;
    for(int i = 0; i<points.size(); i++){
        points[i] = points[i]*center;
    }
}

void Figure3D::rotate(double y, double z) {
    Matrix rotatey;
    rotatey(1,1) = cos(y);
    rotatey(3,3) = cos(y);
    rotatey(1,3) = -sin(y);
    rotatey(3,1) = sin(y);

    Matrix rotatez;
    rotatez(1,1) = cos(z);
    rotatez(2,2) = cos(z);
    rotatez(1,2) = sin(z);
    rotatez(2,1) = -sin(z);

    Matrix m = rotatey*rotatez;

    for(int i = 0; i<points.size(); i++){
        points[i] = points[i]*m;
    }
}

void Figure3D::transform(const double s, double x, double y, double z, Vector3D &v, Vector3D &eye) {
    x = degree_to_rad(x);
    y = degree_to_rad(y);
    z = degree_to_rad(z);
    Matrix scale;
    scale(1,1) = s;
    scale(2,2) = s;
    scale(3,3) = s;

    Matrix rotatex;
    rotatex(2,2) = cos(x);
    rotatex(3,3) = cos(x);
    rotatex(2,3) = sin(x);
    rotatex(3,2) = -sin(x);

    Matrix rotatey;
    rotatey(1,1) = cos(y);
    rotatey(3,3) = cos(y);
    rotatey(1,3) = -sin(y);
    rotatey(3,1) = sin(y);

    Matrix rotatez;
    rotatez(1,1) = cos(z);
    rotatez(2,2) = cos(z);
    rotatez(1,2) = sin(z);
    rotatez(2,1) = -sin(z);

    Matrix center;
    center(4,1) = v.x;
    center(4,2) = v.y;
    center(4,3) = v.z;

    Matrix e;
    double r = sqrt(pow(eye.x,2)+pow(eye.y,2)+pow(eye.z,2));
    double o = atan2(eye.y,eye.x);
    double phi = acos(eye.z/r);

    e(1,1) = -sin(o);
    e(1,2) = -cos(o)*cos(phi);
    e(1,3) = cos(o)*sin(phi);
    e(2,1) = cos(o);
    e(2,2) = -sin(o)*cos(phi);
    e(2,3) = sin(o)*sin(phi);
    e(3,2) = sin(phi);
    e(3,3) = cos(phi);
    e(4,3) = -r;

    Matrix m = rotatex*rotatey*rotatez*center*e;

    for(int i = 0; i<points.size(); i++){
        points[i] = points[i]*scale;
        points[i] = points[i]*m;
    }
}

vector<pair<double,double>> Figure3D::get2DTrianglePoints(double d,double dx, double dy) {
    vector<pair<double,double>> points2D;
    for(int i = 0; i<points.size(); i++) {
        double x = (points[i].x * d) / (-points[i].z);
        double y = (points[i].y * d) / (-points[i].z);
        x+=dx;
        y+=dy;
        points2D.push_back(make_pair(x,y));
    }
    return points2D;
}


vector<Line2D*> Figure3D::get2DLines(bool L3D, double d) {
    img::Color c(ambientReflection[0]*255,ambientReflection[1]*255, ambientReflection[2]*255);
    vector<pair<double,double>> points2D;
    vector<Line2D*> lines;
    for(int i = 0; i<points.size(); i++){
        double x = (points[i].x*d)/(-points[i].z);
        double y = (points[i].y*d)/(-points[i].z);
        points2D.push_back(make_pair(x,y));
        if(first == true){
            xmin = x;
            xmax = x;
            ymin = y;
            ymax = y;
            first = false;
        }
        else{
            if(x<xmin) xmin = x;
            else if(x>xmax) xmax = x;
            if(y<ymin) ymin = y;
            else if(y>ymax) ymax = y;
        }

    }
    for(int i = 0; i<faces.size(); i++){
        if(faces[i].size() == 2){
            Line2D* l = new Line2D();
            l->setPoint1(make_pair(points2D[faces[i][0]].first,points2D[faces[i][0]].second));
            l->setPoint2(make_pair(points2D[faces[i][1]].first,points2D[faces[i][1]].second));
            l->setZ1(points[faces[i][0]].z);
            l->setZ2(points[faces[i][1]].z);
            l->setColor(c);

            lines.emplace_back(l);
        }
        else{
            for(int j = 1; j<faces[i].size(); j++){
                Line2D* l = new Line2D();
                l->setPoint1(make_pair(points2D[faces[i][j-1]].first,points2D[faces[i][j-1]].second));
                l->setPoint2(make_pair(points2D[faces[i][j]].first,points2D[faces[i][j]].second));
                l->setZ1(points[faces[i][j-1]].z);
                l->setZ2(points[faces[i][j]].z);
                l->setColor(c);
                lines.emplace_back(l);
            }
            if(!L3D) {
                Line2D *l = new Line2D();
                l->setPoint1(make_pair(points2D[faces[i][faces[i].size() - 1]].first,points2D[faces[i][faces[i].size() - 1]].second));
                l->setPoint2(make_pair(points2D[faces[i][0]].first, points2D[faces[i][0]].second));
                l->setZ1(points[faces[i][faces[i].size() - 1]].z);
                l->setZ2(points[faces[i][0]].z);
                l->setColor(c);
                lines.emplace_back(l);
            }
        }
    }
    return lines;
}

vector<double> Figure3D::getMinMax() {
    vector<double> minmax;
    minmax.push_back(xmin);
    minmax.push_back(xmax);
    minmax.push_back(ymin);
    minmax.push_back(ymax);
    return minmax;
}

double Figure3D::degree_to_rad(const double degree)
{
    double pi = 3.14159265359;
    return (degree * (pi / 180));
}

vector<vector<int>> Figure3D::getFaces() {
    return faces;
}

Vector3D Figure3D::getPoint(int index) {
    return points[index];
}

int Figure3D::getPointSize() {
    return points.size();
}

Figure3D* Figure3D::makecopy() {
    vector<Vector3D> newpoints = points;
    vector<vector<int>> newfaces = faces;
    Figure3D* f = new Figure3D(newpoints,newfaces, ambientReflection, diffuseReflection,specularReflection,reflectionCoefficient);
    return f;
}

void Figure3D::triangulate() {
    vector<vector<int>> newTrianglefaces;
    for(int i = 0; i<faces.size(); i++){
        for(int j = 1; j<faces[i].size()-1; j++){
            vector<int> newFace = {faces[i][0],faces[i][j],faces[i][j+1]};
            newTrianglefaces.push_back(newFace);
        }
    }
    faces = newTrianglefaces;
}

const vector<double> &Figure3D::getAmbientReflection() const {
    return ambientReflection;
}

const vector<double> &Figure3D::getDiffuseReflection() const {
    return diffuseReflection;
}

const vector<double> &Figure3D::getSpecularReflection() const {
    return specularReflection;
}

double Figure3D::getReflectionCoefficient() const {
    return reflectionCoefficient;
}
