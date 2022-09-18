//
// Created by Niels on 5/03/2021.
//

#ifndef ENGINE_FIGURE3D_H
#define ENGINE_FIGURE3D_H

#include "Line2D.h"
#include "Libraries/vector3d.h"
#include<cmath>
#include<vector>
using namespace std;

class Figure3D {
public:
    Figure3D(const vector<Vector3D> &points, const vector<vector<int>> &faces, const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular = {0.0,0.0,0.0}, const double reflection = 0.0);
    void scale(const double scale);
    void translate(Vector3D &c);
    void rotate(double rotatey, double rotatez);
    void transform(const double scale, double rotatex, double rotatey, double rotatez, Vector3D &center, Vector3D &eye);
    vector<Line2D*> get2DLines(bool b,double d);
    vector<pair<double,double>> get2DTrianglePoints(double d,double dx,double dy);
    vector<double> getMinMax();
    double degree_to_rad(const double degree);
    vector<vector<int>> getFaces();
    Vector3D getPoint(int index);
    int getPointSize();
    Figure3D* makecopy();
    void triangulate();
    const vector<double> &getAmbientReflection() const;
    const vector<double> &getDiffuseReflection() const;
    const vector<double> &getSpecularReflection() const;
    double getReflectionCoefficient() const;
private:
    vector<Vector3D> points;
    vector<vector<int>> faces;
    vector<double> ambientReflection;
    vector<double> diffuseReflection;
    vector<double> specularReflection;
    double reflectionCoefficient;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    bool first = true;

};


#endif //ENGINE_FIGURE3D_H
