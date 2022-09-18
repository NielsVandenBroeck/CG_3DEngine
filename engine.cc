
#define _USE_MATH_DEFINES
#include "Libraries/easy_image.h"
#include "Parsing/ini_configuration.h"
#include "Line2D.h"
#include "Parsing/l_parser.h"
#include "L2D.h"
#include "L3D.h"
#include "Figure3D.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include<cmath>

using namespace std;

img::Figure* makeCube(){
    img::Figure *f =  new img::Figure();
    f->points = {Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),Vector3D::point(-1,-1,1),Vector3D::point(1,1,-1),Vector3D::point(-1,-1,-1),Vector3D::point(1,-1,1),Vector3D::point(-1,1,1)};
    f->faces = {{0,4,2,6},{4,1,7,2},{1,5,3,7},{5,0,6,3},{6,2,7,3},{0,5,1,4}};
    return f;
}

img::Figure* makeTetrahedron(){
    img::Figure *f =  new img::Figure();
    f->points = {Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),Vector3D::point(-1,-1,1)};
    f->faces = {{0,1,2},{1,3,2},{0,3,1},{0,2,3}};
    return f;
}

img::Figure* makeOctahedron(){
    img::Figure *f =  new img::Figure();
    f->points = {Vector3D::point(1,0,0),Vector3D::point(0,1,0),Vector3D::point(-1,0,0),Vector3D::point(0,-1,0),Vector3D::point(0,0,-1),Vector3D::point(0,0,1)};
    f->faces = {{0,1,5},{1,2,5},{2,3,5},{3,0,5},{1,0,4},{2,1,4},{3,2,4},{0,3,4}};
    return f;
}

img::Figure* makeIcosahedron(){
    img::Figure *f =  new img::Figure();
    f->points.reserve(12);
    f->points.push_back(Vector3D::point(0, 0, sqrt(5) / 2));
    for (int i = 0; i < 5; i++) {
        f->points.push_back(Vector3D::point(cos((i * 2 * M_PI) / 5), sin((i * 2 * M_PI) / 5),0.5));
    }
    for (int i = 0; i < 5; i++) {
        f->points.push_back(Vector3D::point(cos((M_PI / 5) + (i * 2 * M_PI / 5)),sin((M_PI / 5) + (i * 2 * M_PI / 5)), -0.5));
    }
    f->points.push_back(Vector3D::point(0, 0, -(sqrt(5) / 2)));
    f->faces = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},{1,6,2},{2,6,7},{2,7,3},{3,7,8},{3,8,4},{4,8,9},{4,9,5},{5,9,10},{5,10,1},{1,10,6},{11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10}};
    return f;
}

img::Figure* makeDodecadron(img::Figure* icosahedron){
    img::Figure *f =  new img::Figure();
    f->points.reserve(20);
    for(int i = 0; i < 20; i++){
        f->points.push_back(Vector3D::point((icosahedron->points[icosahedron->faces[i][0]]+icosahedron->points[icosahedron->faces[i][1]]+icosahedron->points[icosahedron->faces[i][2]])/3));
    }
    f->faces = {{0,1,2,3,4},{0,5,6,7,1},{1,7,8,9,2},{2,9,10,11,3},{3,11,12,13,4},{4,13,14,5,0},{19,18,17,16,15},{19,14,13,12,18},{18,12,11,10,17},{17,10,9,8,16},{16,8,7,6,15},{15,6,5,14,19}};
    return f;
}

img::Figure* makeBuckyBall(img::Figure* ico){
    img::Figure *f =  new img::Figure();
    double mul = 1.0/3.0;
    int count = 0;
    f->points.reserve(120);

    for(auto it = ico->faces.begin(); it != ico->faces.end(); it++){
        f->points.push_back(Vector3D::point((1-mul)*ico->points[(*it)[0]] + mul* ico->points[(*it)[1]]));
        f->points.push_back(Vector3D::point(mul*ico->points[(*it)[0]] + (1-mul)* ico->points[(*it)[1]]));
        f->points.push_back(Vector3D::point((1-mul)*ico->points[(*it)[1]] + mul* ico->points[(*it)[2]]));
        f->points.push_back(Vector3D::point(mul*ico->points[(*it)[1]] + (1-mul)* ico->points[(*it)[2]]));
        f->points.push_back(Vector3D::point((1-mul)*ico->points[(*it)[2]] + mul* ico->points[(*it)[0]]));
        f->points.push_back(Vector3D::point(mul*ico->points[(*it)[2]] + (1-mul)* ico->points[(*it)[0]]));
        f->faces.push_back({count ,count + 1,count + 2,count + 3,count + 4,count + 5});
        count += 6;
    }

    //bovenvlak
    f->faces.push_back({0,6,12,18,24});

    //tussenvlakken
    f->faces.push_back({4,3,33,41,8});
    f->faces.push_back({10,9,48,54,14});
    f->faces.push_back({16,15,57,66,20});
    f->faces.push_back({22,21,69,77,26});
    f->faces.push_back({1,27,81,30,2});

    f->faces.push_back({40,39,91,50,44});
    f->faces.push_back({31,87,94,38,32});
    f->faces.push_back({52,51,97,62,56});
    f->faces.push_back({64,63,103,74,68});
    f->faces.push_back({76,75,109,86,80});

    //ondervlak
    f->faces.push_back({90,96,102,108,114});

    return f;
}

img::Figure* makeCone(const int n, const double h){
    img::Figure *f =  new img::Figure();
    f->points.reserve(n+1);
    f->faces.reserve(n+1);

    vector<int> grondvlak;
    f->points.push_back(Vector3D::point(0,0,h));
    for(int i = 1; i<=n; i++){
        f->points.push_back(Vector3D::point(cos((2*i*M_PI)/n),sin((2*i*M_PI)/n),0));
        grondvlak.push_back(n+1-i);
        f->faces.push_back({i,(i%n)+1, 0});
    }
    f->faces.push_back(grondvlak);
    return f;
}

img::Figure* makeCylinder(const int n, const double h,  bool thicc){
    img::Figure *f =  new img::Figure();
    f->points.reserve(2*n);
    f->faces.reserve(n+2);
    //points
    for(int i = 0; i<n; i++){
        f->points.push_back(Vector3D::point(cos((2*i*M_PI)/n),sin((2*i*M_PI)/n),0));
    }
    for(int i = 0; i<n; i++){
        f->points.push_back(Vector3D::point(cos((2*i*M_PI)/n),sin((2*i*M_PI)/n),h));
    }
    //faces
    for(int i = 0; i<n-1; i++){
        f->faces.push_back({i,(i+1),(n+i+1),n+i});
    }
    f->faces.push_back({n-1,0,n,(2*n)-1});

    if(thicc == false){
        vector<int> grondvlak1;
        vector<int> grondvlak2;
        for(int i = n-1; i>=0; i--){
            grondvlak1.push_back(i);

        }
        for(int i = n; i<2*n; i++){
            grondvlak2.push_back(i);
        }
        f->faces.push_back(grondvlak1);
        f->faces.push_back(grondvlak2);
    }
    return f;
}

img::Figure* makeSphere(const int n){
    img::Figure *f =  makeIcosahedron();
    f->points.reserve(20*pow(3,n)-20);
    f->faces.reserve(20*pow(3,n)-20);
    for(int i = 0; i<n; i++){
        int nrFaces = f->faces.size();
        for(int j = 0; j<nrFaces; j++){
            Vector3D A = f->points[f->faces[j][0]];
            Vector3D B = f->points[f->faces[j][1]];
            Vector3D C = f->points[f->faces[j][2]];
            Vector3D D = Vector3D::point((A.x+B.x)/2 ,(A.y+B.y)/2, (A.z+B.z)/2);
            Vector3D E = Vector3D::point((A.x+C.x)/2 ,(A.y+C.y)/2, (A.z+C.z)/2);
            Vector3D F = Vector3D::point((C.x+B.x)/2 ,(C.y+B.y)/2, (C.z+B.z)/2);
            f->points.push_back(D);
            f->points.push_back(E);
            f->points.push_back(F);

            int pointSize = f->points.size();
            int indexA = f->faces[j][0];
            int indexB = f->faces[j][1];
            int indexC = f->faces[j][2];
            f->faces[j] = {indexA,pointSize-3, pointSize-2};
            f->faces.push_back({indexB, pointSize-1, pointSize-3});
            f->faces.push_back({indexC, pointSize-2, pointSize-1});
            f->faces.push_back({pointSize-3, pointSize-1, pointSize-2});
        }
    }
    for(int i = 0; i<f->points.size(); i++){
        double x = f->points[i].x;
        double y = f->points[i].y;
        double z = f->points[i].z;
        double r = sqrt(pow(x, 2)+pow(y, 2)+pow(z, 2));
        f->points[i] = Vector3D::point(x/r, y/r, z/r);
    }

    return f;
}

img::Figure* makeTorus(const double r, const double R, const int m, const int n){
    img::Figure *f =  new img::Figure();
    f->points.reserve(m*n);
    f->faces.reserve(m*n);
    for(int i=0; i<n; i++){
        double u = (2*i*M_PI)/n;
        for(int j=0; j<m; j++){
            double v = (2*j*M_PI)/m;
            f->points.push_back(Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v)));
            f->faces.push_back({(i*m)+j,(((i+1)%n)*m)+j,((i+1)%n)*m+(j+1)%m, (i*m)+(j+1)%m });
        }
    }
    return f;
}

vector<vector<int>> triangulate(const vector<vector<int>> &faces){
    vector<vector<int>> newTrianglefaces;
    for(int i = 0; i<faces.size(); i++){
        for(int j = 1; j<faces[i].size()-1; j++){
            vector<int> newFace = {faces[i][0],faces[i][j],faces[i][j+1]};
            newTrianglefaces.push_back(newFace);
        }
    }
    return newTrianglefaces;
}

img::EasyImage draw2DLines(const vector<Line2D*> &Lines, int size, const vector<double> &backgroundColor, const vector<double> &minmax, bool zbuffer){
    double xmin = minmax[0];
    double xmax = minmax[1];
    double ymin = minmax[2];
    double ymax = minmax[3];

    //berekend imagex, imagey
    double xrange = xmax-xmin;
    double yrange = ymax-ymin;
    double max;
    if(xrange>=yrange) max = xrange;
    else max = yrange;
    double imagex = size*(xrange/max);
    double imagey = size*(yrange/max);


    //bereken schaalfactor
    double d = 0.95*(imagex/xrange);
    //Middelpunten berekenen
    pair<double,double> MOriginal = make_pair(d*((xmin+xmax)/2), d*((ymin+ymax)/2));
    pair<double, double> MDestination = make_pair(imagex/2,imagey/2);
    double dx = MDestination.first-MOriginal.first;
    double dy = MDestination.second-MOriginal.second;
    //vermenigvuldig alle punten met d
    //punten op juiste positie zetten
    //punten afronden
    for(int i = 0; i < Lines.size(); i++){
        Lines[i]->multiplyPoints(d);
        Lines[i]->addPosition(dx,dy);
        Lines[i]->roundPoints();
    }
    img::Color c(backgroundColor[0]*255,backgroundColor[1]*255,backgroundColor[2]*255);
    img::EasyImage image(round(imagex),round(imagey),c);

    if(zbuffer){
        ZBuffer *zbuff = new ZBuffer(round(imagex), round(imagey));
        for(int i = 0; i < Lines.size(); i++){
            image.draw_zbuff_line(zbuff,round(Lines[i]->getPoint1().first),round(Lines[i]->getPoint1().second),Lines[i]->getZ1(),round(Lines[i]->getPoint2().first),round(Lines[i]->getPoint2().second),Lines[i]->getZ2(),Lines[i]->getColor());
        }
        delete zbuff;
    }
    else{
        for(int i = 0; i < Lines.size(); i++){
            image.draw_line(round(Lines[i]->getPoint1().first),round(Lines[i]->getPoint1().second),round(Lines[i]->getPoint2().first),round(Lines[i]->getPoint2().second),Lines[i]->getColor());
        }
    }
    for(auto it = Lines.begin(); it != Lines.end(); it++){
        delete *it;
    }
    return image;
}

img::EasyImage draw2Dtriangles(vector<Figure3D*> &figures, double size, vector<double> &backgroundColor, const vector<double> &minmax, const vector<Light*> &Lights){
    double xmin = minmax[0];
    double xmax = minmax[1];
    double ymin = minmax[2];
    double ymax = minmax[3];

    //berekend imagex, imagey
    double xrange = xmax-xmin;
    double yrange = ymax-ymin;
    double max;
    if(xrange>=yrange) max = xrange;
    else max = yrange;
    double imagex = size*(xrange/max);
    double imagey = size*(yrange/max);

    //bereken schaalfactor
    double d = 0.95*(imagex/xrange);
    //Middelpunten berekenen
    pair<double,double> MOriginal = make_pair(d*((xmin+xmax)/2), d*((ymin+ymax)/2));
    pair<double, double> MDestination = make_pair(imagex/2,imagey/2);
    double dx = MDestination.first-MOriginal.first;
    double dy = MDestination.second-MOriginal.second;

    img::Color c1(backgroundColor[0]*255,backgroundColor[1]*255,backgroundColor[2]*255);
    img::EasyImage image(round(imagex),round(imagey),c1);
    ZBuffer *zbuff = new ZBuffer(round(imagex), round(imagey));

    vector<Light*> pointlights;
    vector<double> ambientLight = {0.0,0.0,0.0};
    for(int j = 0; j < Lights.size(); j++){
        if(Lights[j]->getType() == "point"){
            pointlights.push_back(Lights[j]);
        }
        ambientLight[0] +=(double) Lights[j]->ambientLight[0];
        ambientLight[1] +=(double) Lights[j]->ambientLight[1];
        ambientLight[2] +=(double) Lights[j]->ambientLight[2];
    }
    if(ambientLight[0] > 1) ambientLight[0] = 1;
    if(ambientLight[1] > 1) ambientLight[1] = 1;
    if(ambientLight[2] > 1) ambientLight[2] = 1;
    if(Lights.empty()){
        ambientLight = {1,1,1};
    }
    for(int i = 0; i < figures.size(); i++){
        Figure3D* Figure = figures[i];
        vector<pair<double,double>> points = Figure->get2DTrianglePoints(d,dx,dy);
        vector<double> ambientReflection = Figure->getAmbientReflection();
        vector<double> diffuseReflection = Figure->getDiffuseReflection();
        vector<double> finalColor1 = {0.0,0.0,0.0};
        finalColor1[0] += ambientLight[0]*ambientReflection[0];
        finalColor1[1] += ambientLight[1]*ambientReflection[1];
        finalColor1[2] += ambientLight[2]*ambientReflection[2];
        vector<vector<int>> faces = Figure->getFaces();
        for(int j = 0; j < faces.size(); j++){
            vector<double> finalColor2 = finalColor1;
            vector<pair<double,double>> newpoints;
            vector<Vector3D> originalpoints;
            for(int k = 0; k < 3; k++){
                newpoints.push_back(points[faces[j][k]]);
                originalpoints.push_back(Figure->getPoint(faces[j][k]));
            }
            if(diffuseReflection[0] != 0.00 || diffuseReflection[1] != 0.00 || diffuseReflection[2] != 0.00){
                Vector3D u = originalpoints[1]-originalpoints[0];
                Vector3D v = originalpoints[2]-originalpoints[0];
                Vector3D w = Vector3D::cross(u,v);
                Vector3D n = Vector3D::normalise(w);
                for(int k = 0; k<Lights.size();k++){
                    if(Lights[k]->getType() == "inf"){
                        vector<double> diffuseLight = Lights[k]->diffuseLight;
                        Vector3D l = Vector3D::normalise(Lights[k]->getlVector());
                        l = -l;
                        double cosa = Vector3D::dot(l,n);
                        if(cosa > 0) {
                            finalColor2[0] += diffuseLight[0] * diffuseReflection[0] * cosa;
                            finalColor2[1] += diffuseLight[1] * diffuseReflection[1] * cosa;
                            finalColor2[2] += diffuseLight[2] * diffuseReflection[2] * cosa;
                        }
                    }
                }
                image.draw_zbuff_triangle(zbuff,newpoints,originalpoints,d, finalColor2, pointlights,dx,dy, diffuseReflection);
            }
            else{
                image.draw_zbuff_triangle(zbuff,newpoints,originalpoints,d, finalColor2);
            }

        }
    }
    for(auto it = figures.begin(); it!= figures.end(); it++){
        delete *it;
    }
    for(auto it = Lights.begin(); it != Lights.end(); it++){
        delete *it;
    }
    delete zbuff;
    return image;
}

void makeMengerSpons(Figure3D* fig, const int nr_iterations, vector<Figure3D*> &figs){
    if(nr_iterations == 0){
        figs.push_back(fig);
    }else{
        vector<Vector3D> locs;
        //alle hoeken krijgen een kubus
        for(int i = 0; i < 8; i++){
            Figure3D* newfig = fig->makecopy();
            newfig->scale(1.0/3.0);
            Vector3D loc = fig->getPoint(i)-newfig->getPoint(i);
            //vector wordt opgeslagen om later posities van andere kubussen te bepalen
            locs.push_back(loc);
            newfig->translate(loc);
            makeMengerSpons(newfig, nr_iterations-1, figs);
        }
        //kubussen tussenin bepalen
        for(int i = 0; i<4; i++){
            for(int j = 4; j < 8; j++){
                if(i+j != 7){
                    Figure3D* newfig = fig->makecopy();
                    newfig->scale(1.0/3.0);
                    Vector3D loc = locs[i]/2+locs[j]/2;
                    newfig->translate(loc);
                    makeMengerSpons(newfig, nr_iterations-1, figs);
                }
            }
        }
        delete fig;
    }
}

void generateFractal(Figure3D* fig, const int nr_iterations, const double scale, vector<Figure3D*> &figs){
    if(nr_iterations == 0){
        figs.push_back(fig);
    }
    else{
        int size = fig->getPointSize();
        for(int i = 0; i < size; i++){
            Figure3D* newfig = fig->makecopy();
            newfig->scale(1/scale);
            Vector3D loc = fig->getPoint(i)-newfig->getPoint(i);
            newfig->translate(loc);
            generateFractal(newfig,nr_iterations-1,scale, figs);
        }
        delete fig;
    }
}

void generateThickFigures(Figure3D* &fig, const double r, const int m, const int n,bool lsystem, vector<Figure3D*> &thickfigures){
    vector<double> ambient = fig->getAmbientReflection();
    vector<double> diffuse = fig->getDiffuseReflection();
    vector<double> specular = fig->getSpecularReflection();

    for(int i = 0; i < fig->getPointSize(); i++){
        img::Figure* bol = makeSphere(m);
        Figure3D* newbol = new Figure3D(bol->points, bol->faces,ambient,diffuse,specular);
        newbol->scale(r);
        Vector3D p = fig->getPoint(i);
        newbol->translate(p);
        thickfigures.push_back(newbol);
    }
    vector<vector<int>> faces = fig->getFaces();
    for(int i = 0; i<faces.size(); i++){
        vector<int> face = faces[i];
        for(int j = 0; j<face.size(); j++){
            if(j == face.size()-1 && lsystem){
                continue;
            }
            Vector3D p1 = fig->getPoint(face[j]);
            Vector3D p2 = fig->getPoint(face[(j+1)%face.size()]);
            if(p1.x == p2.x && p1.y == p2.y && p1.z == p2.z){
                continue;
            }
            Vector3D p1p2 = p2-p1;
            double h = (p1p2.length())/r;
            img::Figure* cilinder = makeCylinder(n,h,true);
            Figure3D* newcilinder = new Figure3D(cilinder->points, cilinder->faces,ambient,diffuse,specular);
            newcilinder->scale(r);
            double R = sqrt(pow(p1p2.x,2)+pow(p1p2.y,2)+pow(p1p2.z,2));
            double o = atan2(p1p2.y,p1p2.x);
            double phi = acos((p1p2.z/R));
            newcilinder->rotate(phi,o);
            Vector3D p = fig->getPoint(face[j]);
            newcilinder->translate(p);
            thickfigures.push_back(newcilinder);
        }
    }
    delete fig;
}

img::EasyImage generate_Image(ini::Configuration &configuration) {
    //2DLSysteem tekenen
    if (configuration["General"]["type"].as_string_or_die() == "2DLSystem"){
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        LParser::LSystem2D l_system;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());
        input_stream >> l_system;
        L2D l = L2D(l_system);
        l.setLineColor(lineColor);
        input_stream.close();
        vector<Line2D *> lines = l.startIteration();
        vector<double> minmax = l.getMinMax();
        return draw2DLines(lines, size, backgroundColor, minmax, false);
    }
    //3D image
    else if(configuration["General"]["type"].as_string_or_die() == "Wireframe" || configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe" || configuration["General"]["type"].as_string_or_die() == "ZBuffering" || configuration["General"]["type"].as_string_or_die() == "LightedZBuffering"){
        bool zBuffer = false;
        bool triangleBuff = false;
        bool Lights = false;
        if(configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){
            zBuffer = true;
        }
        else if(configuration["General"]["type"].as_string_or_die() == "ZBuffering"){
            triangleBuff = true;
        }
        else if(configuration["General"]["type"].as_string_or_die() == "LightedZBuffering"){
            triangleBuff = true;
            Lights = true;
        }
        bool lsystem;
        vector<Line2D*> lines;
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> e = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eye = Vector3D::vector(e[0],e[1],e[2]);
        Matrix meye;
        vector<Light*> lights;
        if(Lights){
            double r = sqrt(pow(eye.x,2)+pow(eye.y,2)+pow(eye.z,2));
            double o = atan2(eye.y,eye.x);
            double phi = acos(eye.z/r);
            meye(1,1) = -sin(o);
            meye(1,2) = -cos(o)*cos(phi);
            meye(1,3) = cos(o)*sin(phi);
            meye(2,1) = cos(o);
            meye(2,2) = -sin(o)*cos(phi);
            meye(2,3) = sin(o)*sin(phi);
            meye(3,2) = sin(phi);
            meye(3,3) = cos(phi);
            meye(4,3) = -r;
            int nrLights = configuration["General"]["nrLights"].as_int_or_die();
            for(int i = 0; i < nrLights; i++){
                vector<double> ambientLight = configuration["Light"+to_string(i)]["ambientLight"].as_double_tuple_or_die();
                vector<double> diffuseLight = configuration["Light"+to_string(i)]["diffuseLight"].as_double_tuple_or_default({0.0,0.0,0.0});
                bool infinity = configuration["Light" + to_string(i)]["infinity"].as_bool_or_default(false);
                if(infinity){
                    vector<double> direction = configuration["Light"+to_string(i)]["direction"].as_double_tuple_or_die();
                    Vector3D d = Vector3D::vector(direction[0],direction[1],direction[2]);
                    d*=meye;
                    infLight* l = new infLight(ambientLight,diffuseLight, Vector3D::vector(d));
                    lights.push_back(l);
                }
                else{
                    if(configuration["Light"+to_string(i)]["location"].exists()){
                        vector<double> location = configuration["Light"+to_string(i)]["location"].as_double_tuple_or_die();
                        double angle = configuration["Light"+to_string(i)]["spotAngle"].as_double_or_default(360);
                        Vector3D loc = Vector3D::point(location[0],location[1],location[2]);
                        loc*=meye;
                        pointLight* l = new pointLight(ambientLight,diffuseLight,loc,angle);
                        lights.push_back(l);
                    }
                    else{
                        Light* l = new Light(ambientLight,diffuseLight);
                        lights.push_back(l);
                    }
                }
            }
        }
        int nFigures = configuration["General"]["nrFigures"].as_int_or_die();
        vector<double> minmax;
        vector<Figure3D*> figures;
        for(int i = 0; i<nFigures; i++){
            double rotatex = configuration["Figure" + to_string(i)]["rotateX"].as_double_or_die();
            double rotatey = configuration["Figure" + to_string(i)]["rotateY"].as_double_or_die();
            double rotatez = configuration["Figure" + to_string(i)]["rotateZ"].as_double_or_die();
            double scale = configuration["Figure" + to_string(i)]["scale"].as_double_or_die();
            vector<double> c = configuration["Figure" + to_string(i)]["center"].as_double_tuple_or_die();
            Vector3D center = Vector3D::vector(c[0],c[1],c[2]);
            vector<double> ambientReflection;
            vector<double> diffuseReflection;
            diffuseReflection = configuration["Figure" + to_string(i)]["diffuseReflection"].as_double_tuple_or_default({0.0,0.0,0.0});
            if(Lights){
                ambientReflection = configuration["Figure" + to_string(i)]["ambientReflection"].as_double_tuple_or_die();
            }
            else{
                ambientReflection = configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die();
            }
            vector<Vector3D> points;
            vector<vector<int>> faces;
            string type = configuration["Figure"+to_string(i)]["type"].as_string_or_die();
            //Lijnen tekenen
            if(type.find("LineDrawing") != string::npos) {
                int FigurePoints = configuration["Figure" + to_string(i)]["nrPoints"].as_int_or_die();
                for (int j = 0; j < FigurePoints; j++) {
                    vector<double> p = configuration["Figure" + to_string(i)]["point" +
                                                                              to_string(j)].as_double_tuple_or_die();
                    Vector3D p3D = Vector3D::point(p[0], p[1], p[2]);
                    points.emplace_back(p3D);
                }
                int FigureFaces = configuration["Figure" + to_string(i)]["nrLines"].as_int_or_die();
                for (int j = 0; j < FigureFaces; j++) {
                    faces.push_back(
                            configuration["Figure" + to_string(i)]["line" + to_string(j)].as_int_tuple_or_die());
                }
            }
                //Kubus tekenen
            else if(type.find("Cube") != string::npos || type == "MengerSponge"){
                img::Figure* f = makeCube();
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Tetrahedron tekenen
            else if(type.find("Tetrahedron") != string::npos){
                img::Figure* f = makeTetrahedron();
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Octahedron tekenen
            else if(type.find("Octahedron") != string::npos){
                img::Figure* f = makeOctahedron();
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Icosahedron tekenen
            else if(type.find("Icosahedron") != string::npos) {
                img::Figure* f = makeIcosahedron();
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Dodecahedron tekenen
            else if(type.find("Dodecahedron") != string::npos){
                img::Figure* f = makeDodecadron(makeIcosahedron());
                points = f->points;
                faces = f->faces;
                delete f;
            }
            else if(type.find("BuckyBall") != string::npos){
                img::Figure* f = makeBuckyBall(makeIcosahedron());
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Cone tekenen
            else if(type =="Cone"){
                int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
                img::Figure* f = makeCone(configuration["Figure"+to_string(i)]["n"].as_int_or_die(), configuration["Figure"+to_string(i)]["height"].as_double_or_die());
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Cylinder tekenen
            else if(type =="Cylinder"){
                int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
                img::Figure* f = makeCylinder(configuration["Figure"+to_string(i)]["n"].as_int_or_die(), configuration["Figure"+to_string(i)]["height"].as_double_or_die(),false);
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Sphere tekenen
            else if(type =="Sphere"){
                int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
                img::Figure* f = makeSphere(configuration["Figure"+to_string(i)]["n"].as_int_or_die());
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //Torus tekenen
            else if(type =="Torus"){
                int n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
                img::Figure* f = makeTorus(configuration["Figure"+to_string(i)]["r"].as_double_or_die(),configuration["Figure"+to_string(i)]["R"].as_double_or_die(),configuration["Figure"+to_string(i)]["m"].as_int_or_die(),configuration["Figure"+to_string(i)]["n"].as_int_or_die());
                points = f->points;
                faces = f->faces;
                delete f;
            }
                //3D L-systeem tekenen
            else if(type.find("3DLSystem") != string::npos){
                LParser::LSystem3D l_system;
                std::ifstream input_stream(configuration["Figure"+to_string(i)]["inputfile"].as_string_or_die());
                input_stream >> l_system;
                L3D l = L3D(l_system);
                input_stream.close();
                img::Figure* f = l.startIteration();
                points = f->points;
                faces = f->faces;
                lsystem = true;
                delete f;
            }
            if(triangleBuff && type.find("Thick") == string::npos){
                faces = triangulate(faces);
            }
            Figure3D* f = new Figure3D(points, faces,ambientReflection,diffuseReflection);
            if(type.find("Thick") != string::npos){
                double r = configuration["Figure"+to_string(i)]["radius"].as_double_or_die();
                double m = configuration["Figure"+to_string(i)]["m"].as_int_or_die();
                double n = configuration["Figure"+to_string(i)]["n"].as_int_or_die();
                vector<Figure3D*> thickFigures = {};
                generateThickFigures(f, r, m, n,lsystem, thickFigures);
                for(int j = 0; j < thickFigures.size(); j++){
                    if(triangleBuff){
                        thickFigures[j]->triangulate();
                    }
                    thickFigures[j]->transform(scale, rotatex, rotatey, rotatez, center, eye);
                    vector<Line2D *> lines1 = thickFigures[j]->get2DLines(lsystem,1);
                    lines.insert(lines.end(), lines1.begin(), lines1.end());
                    if (i == 0 && j == 0) {
                        minmax = thickFigures[j]->getMinMax();
                    } else {
                        vector<double> newminmax = thickFigures[j]->getMinMax();
                        if (newminmax[0] < minmax[0]) minmax[0] = newminmax[0];
                        if (newminmax[1] > minmax[1]) minmax[1] = newminmax[1];
                        if (newminmax[2] < minmax[2]) minmax[2] = newminmax[2];
                        if (newminmax[3] > minmax[3]) minmax[3] = newminmax[3];
                    }
                }
                figures.insert(figures.end(),thickFigures.begin(),thickFigures.end());
            }
            else if(type == "MengerSponge"){
                vector<Figure3D*> figs = {};
                makeMengerSpons(f,configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die(), figs);
                for(int j = 0; j < figs.size(); j++){
                    figs[j]->transform(scale, rotatex, rotatey, rotatez, center, eye);
                    vector<Line2D *> lines1 = figs[j]->get2DLines(lsystem,1);
                    lines.insert(lines.end(), lines1.begin(), lines1.end());
                    if (i == 0 && j ==0) {
                        minmax = figs[j]->getMinMax();
                    } else {
                        vector<double> newminmax = figs[j]->getMinMax();
                        if (newminmax[0] < minmax[0]) minmax[0] = newminmax[0];
                        if (newminmax[1] > minmax[1]) minmax[1] = newminmax[1];
                        if (newminmax[2] < minmax[2]) minmax[2] = newminmax[2];
                        if (newminmax[3] > minmax[3]) minmax[3] = newminmax[3];
                    }
                }
                figures.insert(figures.end(),figs.begin(),figs.end());
            }
            else if(type.find("Fractal") != string::npos){
                vector<Figure3D*> fractFigures = {};
                generateFractal(f,configuration["Figure"+to_string(i)]["nrIterations"].as_int_or_die(),configuration["Figure"+to_string(i)]["fractalScale"].as_double_or_die(), fractFigures);
                for(int j = 0; j < fractFigures.size(); j++){
                    fractFigures[j]->transform(scale, rotatex, rotatey, rotatez, center, eye);
                    vector<Line2D *> lines1 = fractFigures[j]->get2DLines(lsystem,1);
                    lines.insert(lines.end(), lines1.begin(), lines1.end());
                    if (i == 0 && j ==0) {
                        minmax = fractFigures[j]->getMinMax();
                    } else {
                        vector<double> newminmax = fractFigures[j]->getMinMax();
                        if (newminmax[0] < minmax[0]) minmax[0] = newminmax[0];
                        if (newminmax[1] > minmax[1]) minmax[1] = newminmax[1];
                        if (newminmax[2] < minmax[2]) minmax[2] = newminmax[2];
                        if (newminmax[3] > minmax[3]) minmax[3] = newminmax[3];
                    }
                }
                figures.insert(figures.end(),fractFigures.begin(),fractFigures.end());
            }
            else{
                f->transform(scale, rotatex, rotatey, rotatez, center, eye);
                figures.emplace_back(f);
                vector<Line2D *> lines1 = f->get2DLines(lsystem,1);
                lines.insert(lines.end(), lines1.begin(), lines1.end());
                if (i == 0) {
                    minmax = f->getMinMax();
                } else {
                    vector<double> newminmax = f->getMinMax();
                    if (newminmax[0] < minmax[0]) minmax[0] = newminmax[0];
                    if (newminmax[1] > minmax[1]) minmax[1] = newminmax[1];
                    if (newminmax[2] < minmax[2]) minmax[2] = newminmax[2];
                    if (newminmax[3] > minmax[3]) minmax[3] = newminmax[3];
                }
            }
        }
        if(triangleBuff){
            for(auto it = lines.begin(); it != lines.end(); it++){
                delete *it;
            }
            return draw2Dtriangles(figures, size, backgroundColor, minmax, lights);
        }
        else{
            for(auto it = figures.begin(); it!= figures.end(); it++){
                delete *it;
            }
            return draw2DLines(lines,size,backgroundColor,minmax,zBuffer);
        }
    }
    else{
        cerr<<"this type has not been implemented."<<endl;
        return img::EasyImage(0,0);
    }

}

int main(int argc, char const* argv[])
{
    auto start = std::chrono::system_clock::now();
    int retVal = 0;
    try
    {
        std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string fileName;
            while (std::getline(fileIn, fileName)) {
                args.push_back(fileName);
            }
        }
        for(std::string fileName : args)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }
            img::EasyImage image = generate_Image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                //std::string fileName(fileName);
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;
                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    return retVal;
}