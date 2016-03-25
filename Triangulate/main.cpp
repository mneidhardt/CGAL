//
//  main.cpp
//  Triangulate point set using CGAL.
//
//  Created by michael neidhardt on 05/04/15.
//  Copyright (c) 2015 michael neidhardt. All rights reserved.
//

#include <iostream>
#include <iomanip>   // for flush
#include <cmath>
#include <GLUT/glut.h>
#include <vector>
#include <fstream>
#include <string>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include "Datasource2.h"

#include <fstream>
#include <cassert>
// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
// typedefs for defining the adaptor
typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K2;
typedef CGAL::Delaunay_triangulation_2<K2>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;
// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;
typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle2;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef VD::Edge_iterator             EdgeIterator;



typedef double                                                 Coord_type;
typedef CGAL::Simple_cartesian<Coord_type>                     K;
typedef K::Point_2                                             dtPoint;
typedef CGAL::Triangulation_vertex_base_2<K>                   Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>       Vb;
typedef CGAL::Triangulation_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>            Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                  Dtr;
typedef CGAL::Triangulation_hierarchy_2<Dtr>                   Delaunay;
/* ---------------------------------------------------------------- */

typedef Delaunay::Vertex_iterator           Vertex_iterator;
typedef Delaunay::Vertex_circulator         Vertex_circulator;
typedef Delaunay::Vertex_handle             Vertex_handle;

using namespace std;

typedef float T;

/*double fillItUp(Delaunay& dt, int argc, char** argv);
 void findAllNearestNeighbors(Delaunay& dt, vector<Point>& pointset,
 vector<pair<Point, Point> >& ann);
 */

VD vd;
Delaunay dt;
vector<dtPoint> dpoints;
vector<Site_2>  vdpoints;

float translate_x = 560000.0f;
float translate_y = 6129000.0f;
float translate_z = 0.0f;
float rotangle = 0.0f;
float eyeX = 0.0f;
float eyeY = 1.5f;
float eyeZ = 0.0f;
float lookX = 0.0f;
float lookY = 0.0f;
float lookZ = 0.0f;
GLfloat pointsize = 4.0;
GLint rectsize=700;
int p = 1;  // What to show.


void init(void) {
    GLfloat values[2];
    glGetFloatv (GL_LINE_WIDTH_GRANULARITY, values);
    glGetFloatv (GL_LINE_WIDTH_RANGE, values);
    glEnable (GL_LINE_SMOOTH);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glLineWidth (0.3);
    
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glColor3f(0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
}

void drawLine(Point_2 p, Point_2 q) {
    glBegin(GL_LINES);
    glVertex3f(p.x(), p.y(), 0.0);
    glVertex3f(q.x(), q.y(), 0.0);
    glEnd();
}

void drawPoint(Point_2 p) {
    glPointSize(pointsize);
    glBegin(GL_POINTS);
    glVertex3f(p.x(), p.y(), 0.0);
    glEnd();
}

void print_endpoint(Halfedge_handle e, bool is_src) {
    std::cout << "\t";
    if ( is_src ) {
        if ( e->has_source() ) {
            std::cout << e->source()->point() << std::endl;
            drawPoint(e->source()->point());
        }
        else {
            std::cout << "point at infinity" << std::endl;
        }
    } else {
        if ( e->has_target() ) {
            std::cout << e->target()->point() << std::endl;
            drawPoint(e->source()->point());
        }
        else {
            std::cout << "point at infinity" << std::endl;
        }
    }
}

void pointLocation(Point_2 p) {
    std::cout << "Query point (" << p.x() << "," << p.y()
    << ") lies on a Voronoi " << std::flush;
    Locate_result lr = vd.locate(p);
    if ( Vertex_handle* v = boost::get<Vertex_handle>(&lr) ) {
        std::cout << "vertex." << std::endl;
        std::cout << "The Voronoi vertex is:" << std::endl;
        std::cout << "\t" << (*v)->point() << std::endl;
    } else if ( Halfedge_handle* e = boost::get<Halfedge_handle>(&lr) ) {
        std::cout << "edge." << std::endl;
        std::cout << "The source and target vertices "
        << "of the Voronoi edge are:" << std::endl;
        print_endpoint(*e, true);
        print_endpoint(*e, false);
    } else if ( Face_handle* f = boost::get<Face_handle>(&lr) ) {
        std::cout << "face." << std::endl;
        std::cout << "The vertices of the Voronoi face are"
        << " (in counterclockwise order):" << std::endl;
        Ccb_halfedge_circulator ec_start = (*f)->ccb();
        Ccb_halfedge_circulator ec = ec_start;
        
        Point_2 prev = ec->source()->point();
        ++ec;
        while (ec != ec_start) {
            glBegin(GL_LINES);
            glVertex3f(prev.x(), prev.y(), 0.0);
            glVertex3f(ec->source()->point().x(), ec->source()->point().y(), 0.0);
            glEnd();
            prev = ec->source()->point();
            ++ec;
        }
    }
    std::cout << std::endl;
}

void drawAllPoints() {
    for (int i=0; i<vdpoints.size(); i++) {
        drawPoint(vdpoints[i]);
    }
}

void drawTriangulation() {
    Vertex_iterator vi;
    
    // Outer loop runs over all nodes of the triangulation.
    for (vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) {
        Vertex_circulator vc = dt.incident_vertices(vi);
        if (dt.is_infinite(vc)) { ++vc; }
        dtPoint firstIncidentPoint(vc->point());
        
        // Inner loop runs over all edges for each node.
        do {
            if ( ! dt.is_infinite(vc)) {   //OR THIS?: if (vc != dt.infinite_vertex())
                glBegin(GL_LINES);
                glVertex3f(vi->point().x(), vi->point().y(), 0.0);
                glVertex3f(vc->point().x(), vc->point().y(), 0.0);
                glEnd();
            }
            ++vc;
        } while (vc->point() != firstIncidentPoint);
    }
}

void drawVoronoi() {
    for (EdgeIterator eit = vd.edges_begin(); eit!=vd.edges_end(); ++eit) {
        //drawLine(eit->source()->point(), eit->target()->point());
        if (eit->has_source() && eit->has_target()) {
            drawLine(eit->source()->point(), eit->target()->point());
        }
    }
    
/*    for (VD::Face_iterator fit=vd.faces_begin(),fit_end=vd.faces_end();
         fit!=fit_end;
         ++fit) {
        
        Ccb_halfedge_circulator ccb = fit->ccb();
        // Nu skal jeg så finde ud af at iterere over denne circulator
        // og for hver halfedge, tegne den linie den repræsenterer.
        glVertex3f(fit->dual()->point().x(), fit->dual()->point().y(), 0.0);
    }
    glEnd();
*/
    /*for (int i=0; i<vdpoints.size(); i++) {
        pointLocation(vdpoints[i]);
    }*/
    
}


void display(void) {
    glClear (GL_COLOR_BUFFER_BIT);
    
    glPushMatrix();

    if ((p&1) > 0) {
        drawAllPoints();
    }

    if ((p&2) > 0) {
        drawTriangulation();
    }
    
    if ((p&4) > 0) {
        drawVoronoi();
    }
    glPopMatrix();
    
    glutSwapBuffers();
}



void reshape (int w, int h) {
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    
    if (w <= h)
        gluOrtho2D(0.0, GLfloat(rectsize),
                   0.0, GLfloat(rectsize)*(GLfloat) h/(GLfloat) w);
    else
        gluOrtho2D(0.0, GLfloat(rectsize)*(GLfloat) w/(GLfloat) h,
                   0.0, GLfloat(rectsize));
    
    glMatrixMode(GL_MODELVIEW);
    
}


void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:
            exit(0);
            break;
    }
}

int main(int argc, char** argv) {
    
    if (argc < 2) {
        std::cout << "Syntaxen er: " << argv[0] << " n r [p]\n" <<
                  "n = number of points\n" <<
                  "r = window size\n" <<
                  "p: 1 = draw points\n" <<
                  "   2 = draw Delaunay\n" <<
                  "   3 = draw Delaunay and points\n" <<
                  "   4 = draw Voronoi\n" <<
                  "   5 = draw Voronoi and points\n" <<
                  "   6 = draw Voronoi and Delaunay\n" <<
                  "   7 = draw Voronoi, Delaunay and points\n";
        
        return 0;
    }

    if (argc > 2) { rectsize = atoi(argv[2]); }
    if (argc > 3) { p = atoi(argv[3]); }

    if ((p&1) > 0) {
        cout << "points ";
    }
    if ((p&2) > 0) cout << "Delaunay ";
    if ((p&4) > 0) cout << "Voronoi";
    cout << "\n";
    
    int n = atoi(argv[1]);
    Datasource2<Coord_type, dtPoint> dsrc;
    dpoints = dsrc.getRandomDoublePoints(n, true, rectsize);
    
    // Transform dpoints to a pointset that VD accepts:
    for (int i=0; i<dpoints.size(); i++) {
        Site_2 t(dpoints[i].x(), dpoints[i].y());
        vdpoints.push_back(t);
    }
    
    dt.insert(dpoints.begin(), dpoints.end());
    vd.insert(vdpoints.begin(), vdpoints.end());
    
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(rectsize, rectsize);
    glutInitWindowPosition (100, 100);
    glutCreateWindow (argv[0]);
    init ();
    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutReshapeFunc(reshape);
    //glutMotionFunc(mousemotionfunc);
    //glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}

