/* This is a library that can create sets of points for use in testing my ANN
+  programs.
+  It is based on Datasource.h, which was getting messy so I moved it to this.
+  I need points with different distributions and types.
+  It is all 2D. This is supposed to be used to create files of point sets.
+
+  The class is parametrized on coordinate type and node (point) type.
+  I.e. you must supply the type for the individual coordinate, T, and
+  the type for the node (point), N.
+  With regards to rational numbers, this only supports leda::rational directly.
+  Based on Pointsource.h, and made initially for the ANN problem.
+
+  Michael Neidhardt, august 2008.
*/

#ifndef DATASOURCE2_H
#define DATASOURCE2_H

#include <fstream>
#include <string>
// I do not have Leda anymore...
// #include <LEDA/core/random.h>
// #include <LEDA/numbers/integer.h>
#include <vector>
#include <math.h>

#define PI 3.141592926535897932384626433832

using namespace std;

template<typename T, typename N>
class Datasource2 {
public:
   Datasource2();

   /* These return vectors of points/pointers to vectors of points: */
   /* First the floating point generators: */
   vector<N>   getRandomDoublePoints(unsigned int n, bool onUnitSquare=false, float factor=1.0);
   vector<N>   getRDInCircle(unsigned int n, float radius=1.0);
   vector<N>   getRDOnUnitCircle(unsigned int n);
   vector<N>   getRDOnLine(unsigned int n);
   /* Here the integral point generators: */
   // vector<N>   getRandomIntPoints(unsigned int n, int precision=31);
   vector<N>   getRandomABDPYPoints(unsigned int n);
   vector<N>   getRIInUnitCircle(unsigned int n);
   vector<N>   getRIOnUnitCircle(unsigned int n);

   vector<N>   getKMPSY1(bool includebig=true);
   vector<N*>* getTTLKMPSY1(bool includebig=true);

   /* These return either a single point or a single number: */
   double doublerand();
   double uniformdoublerand();
   double uniformdoublerandUNITSQ();
   double randomInt(int exponent); 
   N      randomPointInCircle(T radius=1, bool keepdecimals=true);
   N      randomPointOnCircle(T radius=1, bool keepdecimals=true);

   /* File i/o related functions. */
   void        writePointset(vector<N>& pointset, const char* outfile);
   void        writePointset4Triangle(vector<N>& pointset, const char* outfile);
   // void        writeLedaRationalPointset(vector<N>& pointset, const char* outfile);
   vector<N>   readPointset(const char* infile);
   vector<N>   readPointset(const char* infile, const int i);
   vector<N>   readLedaRationalPointset(const char* infile);
   vector<N*>* readTTLPointset(const char* infile);

   /* Some utility functions... */
   void print(vector<N>& nodes);
   void printPS(vector<N>& pointset);
    void translatePointset(vector<N>& pointset, float translX, float translY);
    void dumpPointset(vector<N>& pointset);

};

// -------- Implementation here: ----------

template<typename T, typename N>
Datasource2<T,N>::Datasource2() { srand(1); }


/* Returns n random points with double coordinates. Not used so far - 
+  might cause overflow as the exponent can be up to 1023.
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRandomDoublePoints(unsigned int n, bool onUnitSquare, float factor) {

   vector<N> pointset;

   for (unsigned int i=0; i<n; i++) {
      if (onUnitSquare) {
         pointset.push_back(N(factor*uniformdoublerandUNITSQ(), factor*uniformdoublerandUNITSQ()));
      }
      else {
         pointset.push_back(N(uniformdoublerand(), uniformdoublerand()));
      }
   }

   return pointset;
}


/* Returns n random double points inside circle with given radius. Default radius is 1. */
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRDInCircle(unsigned int n, float radius) {
    
    vector<N> pointset;
    
    for (unsigned int i=0; i<n; i++) {
        pointset.push_back(randomPointInCircle(radius));
    }
    
    return pointset;
}



/* Returns n random double points approx. on the unit circle. 
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRDOnUnitCircle(unsigned int n) {

   vector<N> pointset;

   for (unsigned int i=0; i<n; i++) {
      pointset.push_back(randomPointOnCircle(1));
   }

   return pointset;
}

/* Returns n random double points approx. on a horizontal line. 
   This was made to be nasty to KDTrees.
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRDOnLine(unsigned int n) {

   vector<N> pointset;

   for (unsigned int i=0; i<n; i++) {
      pointset.push_back(
               N(uniformdoublerand(), uniformdoublerandUNITSQ())
               );
   }

   return pointset;
}


// Returns n random points with integer coordinates:
/*
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRandomIntPoints(unsigned int n, int precision) {

   vector<N> pointset(0);

   leda::random_source S(precision);
   int x,y;
   
   for (unsigned int i=0; i<n; i++) {
      S >> x >> y;
      pointset.push_back(N(T(x), T(y)));
   }

   return pointset;

}
*/

/* Returns n random points with integral coords, but stored in doubles.
+  The coords are in the range ]-2^51, 2^51[.
+  This is for use with TTL and ABDPY traits.
*/
template <typename T, typename N>
vector<N> Datasource2<T,N>::getRandomABDPYPoints(unsigned int n) {
    
   vector<N> pointset(0);

   for (unsigned int i=0; i<n; ++i) {
      pointset.push_back(N(randomInt(51), randomInt(51)));
   }

   return pointset;
}


/* Returns n random integer points inside a circle with radius 1000. */
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRIInUnitCircle(unsigned int n) {

   vector<N> pointset;

   for (unsigned int i=0; i<n; i++) {
      pointset.push_back(randomPointInCircle(1000, false));
   }

   return pointset;
}

/* Returns n random integer points approx. on a circle with radius 1000. */
template<typename T, typename N>
vector<N> Datasource2<T,N>::getRIOnUnitCircle(unsigned int n) {

   vector<N> pointset;

   for (unsigned int i=0; i<n; i++) {
      pointset.push_back(randomPointOnCircle(1000, false));
   }

   return pointset;
}

/* Returns 256*128 points that are nearly collinear. The idea is from the article "Classroom
+  examples of robustness problems in geometric computations" by
+  Kettner,Mehlhorn,Pion,Schirra,Yap - 2006.
+  This function returns data with 2 fixed points, and one other that is varied
+  similar to the one that generates the picture in figure 2(a) on page 5. 
+  The argument decides whether to include the two 'big numbered' points or not.
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::getKMPSY1(bool includebig) {

   vector<N> pointset;

   if (includebig) {
      pointset.push_back(N(12.0, 12.0));
      pointset.push_back(N(24.0, 24.0));
   }

   double eps = pow(2.0, -53);

   for (unsigned int x=0; x<256; x++) {
      for (unsigned int y=0; y<=x; y++) {
         pointset.push_back(N(0.5+x*eps, 0.5+y*eps));
      }
   }

/*
   std::cout.precision(53);
   std::cout.setf(ios::fixed);
   for (unsigned int i=0; i<pointset.size(); i++) {
      std::cout << i << " " << pointset[i].x() << ", " << pointset[i].y() << "\n";
   }
*/

   return pointset;
}


/* Since TTL expects a pointer to a vector of points, this is what I 
+  create here. All it does is repackage a vector into a pointer to a vector of pointers.
*/
template<typename T, typename N>
vector<N*>* Datasource2<T,N>::getTTLKMPSY1(bool includebig) {

   vector<N> pointset = getKMPSY1(includebig);
   unsigned int pidx=0;
   T x, y;

   vector<N*>* points = new vector<N*>(pointset.size());
   typename vector<N*>::iterator it;

   for (it = points->begin(); it != points->end(); ++it) {
      T x = pointset[pidx].x();
      T y = pointset[pidx].y();
      *it = new N(x, y);
      ++pidx;
   }

   return points;
}



/**** Below are functions that create a random number or a random point in 2D. */

/*****************************************************************************/
/*                                                                           */
/*  doublerand()   Generate a double with random 53-bit significand and a    */
/*                 random exponent in [0, 511]. By J.R. Shewchuk             */
/*                                                                           */
/* NB: I have decided this doesnt work for me - dont know exactly what's     */
/* wrong, but even triangle segfaults when given data from this (saved    */
/* to file by me). I use uniformdoublerand in stead.                         */
/*****************************************************************************/
template<typename T, typename N>
double Datasource2<T,N>::doublerand() {
   double result;
   double expo;
   long a, b, c;
   long i;

   a = random();
   b = random();
   c = random();
   result = (double) (a - 1073741824) * 8388608.0 + (double) (b >> 8);
   for (i = 512, expo = 2; i <= 131072; i *= 2, expo = expo * expo) {
      if (c & i) {
         result *= expo;
      }
   }
   return result;
}

/*****************************************************************************
+                                                                            
+  uniformdoublerand()   Generate a double with random 53-bit significand
+  in the range [-2^53, 2^53-1].
+ By J.R.Shewchuk.            
***************************************************************************** */
template<typename T, typename N>
double Datasource2<T,N>::uniformdoublerand() {
   double result;
   long a, b;

   a = random();
   b = random();
   result = (double) (a - 1073741824) * 8388608.0 + (double) (b >> 8);
   
   return result;
}

/******************************************************************************/
/*                                                                            */
/*  uniformdoublerandUNSQ() Generate a random double on the unit square.      */
/*                                                                            */
/******************************************************************************/
template<typename T, typename N>
double Datasource2<T,N>::uniformdoublerandUNITSQ() {
   /*
   double result;
   long a, b;

   a = random();
   if (a < 1073741824) { a = 1073741824; }
   b = random();
   result = (double) (a - 1073741824) * 8388608.0 + (double) (b >> 8);
   result = result/((RAND_MAX - 1073741824) * 8388608.0 + (double) (RAND_MAX >> 8));
   */

   // The above gives to high a density along the axes.
   // The following works, but does not have that many decimals, though.
   return double(double(random())/double(RAND_MAX));
}


/* Returns an integer stored in a double, in the interval ]-2^exponent, 2^exponent[.
*/
template<typename T, typename N>
double Datasource2<T,N>::randomInt(int exponent) {
   double unitrand = double(random())/double(RAND_MAX);
   
   double x = unitrand*(-1.0 + pow(2.0, exponent+1)) - (-1.0 + pow(2.0, exponent));

   return floor(x);
}



/* Returns a random point inside a circle of the given radius, 
+  which is 1 by default. See text on the method at the bottom of this file.
+  If keepdecimals is false (the default is true), the coords are rounded
+  in order to drop the decimals. This so it can be used by programs expecting
+  integer input. It is however a crude method, and does not look very random!!!
+  If keepdecimals, the coordinates keep the decimals.
*/
template<typename T, typename N>
N Datasource2<T,N>::randomPointInCircle(T radius, bool keepdecimals) {
   // First get a uniformly random number in interval [0, 2PI].
   double theta = (2*PI)*double(random())/double(RAND_MAX);
   // Then get 2 numbers, both uniformly random in [0, radius]:
   double r1 = double(radius)*double(random())/double(RAND_MAX);
   double r2 = double(radius)*double(random())/double(RAND_MAX);
   // Then get the biggest of those 2, producing M, the length of the vector.
   // This should give me a 'ramp' distribution along the vector, i.e. I avoid
   // too high a density towards the center:
   double M = max(r1, r2);

   if (keepdecimals) { return N(M*cos(theta), M*sin(theta)); }
   else              { return N(round(M*cos(theta)), round(M*sin(theta))); }
}

/* Returns a random point approximately(!) on a circle of the given radius, 
+  which is 1 by default. See text at bottom of file.
+  This creates points with double as coordinates.
*/
template<typename T, typename N>
N Datasource2<T,N>::randomPointOnCircle(T radius, bool keepdecimals) {
   // First get a uniformly random number in interval [0, 2PI].
   double theta = (2*PI)*double(random())/double(RAND_MAX);

   // spread is how much on either side of the circle the point can drift.
   // The bigger the spread, the smaller the 'band' the point will fall in.
   // Dont mind a bit of round off error here, so dont use a power of 2.
   double spread=1000.0;
   // Create M as random number approx. in [radius-radius/spread, radius+radius/spread].
   double M = (2.0*radius/spread)*double(random())/RAND_MAX;
   M = M+(radius-radius/spread);

   if (keepdecimals) { return N(M*cos(theta), M*sin(theta)); }
   else              { return N(round(M*cos(theta)), round(M*sin(theta))); }
}


/**** End of random number/point functions. ***/


/********************************************************************************/
/**** Below are functions for reading and writing point sets. ****/

/* Write a pointset to a file. */
template<typename T, typename N>
void Datasource2<T,N>::writePointset(vector<N>& pointset, const char* outfile) {
   ofstream outStream(outfile);
   outStream.precision(53);
   outStream.setf(ios::fixed);
   //outStream.setf(ios::showpoint);

   if (outStream.fail()) {
      cerr << "Cannot open file: ";
      printf("%s\n", outfile);
      exit(1);
   }

   outStream << pointset.size() << "\n";
   
   for (unsigned int i=0; i<pointset.size(); i++) {
      outStream << pointset[i].x() << " " << pointset[i].y() << "\n";
   }

   outStream.close();
}


/* Write a pointset to a file, to be used by Triangle by Shewchuk.
   Almost similar to writePointset. */
template<typename T, typename N>
void Datasource2<T,N>::writePointset4Triangle(vector<N>& pointset, const char* outfile) {
   ofstream outStream(outfile);
   outStream.precision(53);
   outStream.setf(ios::fixed);
   //outStream.setf(ios::showpoint);

   if (outStream.fail()) {
      cerr << "Cannot open file: ";
      printf("%s\n", outfile);
      exit(1);
   }

   outStream << pointset.size() << " 2 0 0\n";
   
   for (unsigned int i=0; i<pointset.size(); i++) {
      outStream << (i+1) << " " << pointset[i].x() << " " << pointset[i].y() << "\n";
   }

   outStream.close();
}

/* This assumes that T equals leda::rational, or that T has a member function numerator().
+  It writes only the numerator, as it assumes that the denominator is 1. This was originally
+  made for use with integer input, and so this is an easy way to store rationals. */
/*
template<typename T, typename N>
void Datasource2<T,N>::writeLedaRationalPointset(vector<N>& pointset, const char* outfile) {
   ofstream outStream(outfile);
   outStream.precision(53);
   outStream.setf(ios::fixed);
   outStream << pointset.size() << "\n";
   
   for (unsigned int i=0; i<pointset.size(); i++) {
      outStream << pointset[i].x().numerator() << " " 
                << pointset[i].y().numerator() << "\n";
   }

   outStream.close();
}
*/


/* Reads a pointset from a file into a point/node with coordinate type T.
   I expect a file with 3D points, with x,y,z coordinates on each line, separated by space,
   like so:
   1.2 3.3 2.0
   2.7 8.7 4.1
   ...
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::readPointset(const char* infile) {
   return readPointset(infile, 1);
}

/* Reads a pointset from a file into a point/node with coordinate type T.
   I expect a file with 3D points, with x,y,z coordinates on each line, separated by space,
   like so:
   1.2 3.3 2.0
   2.7 8.7 4.1
   ...

   The second arg says how many to skip - if i is 2, only read in every other line,
   and if i is 3, read in every 3rd line.
*/
template<typename T, typename N>
vector<N> Datasource2<T,N>::readPointset(const char* infile, const int i) {
   vector<N> pointset(0);
   
   ifstream inStream(infile);
   inStream.precision(53);
   inStream.setf(ios::fixed);
   
   if (inStream.fail()) {
      cerr << "Cannot open file.\n";
      exit(1);
   }

   T x, y, z;
   long counter = 0;

   while (inStream >> x >> y >> z) {
       ++counter;
       if (counter%i == 0) {
           pointset.push_back(N(x, y, z));
       }
   }

   inStream.close();
   
   return pointset;
}

/* When reading data into leda::rationals, I assume they are stored as only numerators,
+  so for 2D points (one per line) I expect 2 numbers per line. These can be of type
+  int or double.
+  They get turned into leda::rationals which in turn get made into a point, that
+  gets appended to pointset. */
template<typename T, typename N>
vector<N> Datasource2<T,N>::readLedaRationalPointset(const char* infile) {
   vector<N> pointset(0);
   
   ifstream inStream(infile);
   inStream.precision(53);
   inStream.setf(ios::fixed);
   
   if (inStream.fail()) {
      cerr << "Cannot open file.\n";
      exit(1);
   }

   int numOfPoints;
   double x, y;
   inStream >> numOfPoints;     // Used to use the number of points, which is
                                // is on first line in file, but with
                                // push_back it is not necessary...

   while (inStream >> x >> y) {
      pointset.push_back(N(T(x), T(y)));
   }

   inStream.close();
   
   return pointset;
}

/* Since TTL expects a pointer to a vector of points, this is what I do
+  here. */
template<typename T, typename N>
vector<N*>* Datasource2<T,N>::readTTLPointset(const char* infile) {

   ifstream inStream(infile);
   inStream.precision(53);
   inStream.setf(ios::fixed);
      
   if (inStream.fail()) {
      cerr << "Cannot open file.\n";
      exit(1);
   }
      
   int numOfPoints;
   T x, y;
   inStream >> numOfPoints;

   vector<N*>* points = new vector<N*>(numOfPoints);
   typename vector<N*>::iterator it;

   for (it = points->begin(); it != points->end(); ++it) {
      inStream >> x >> y;
      *it = new N(x, y);
   }

   inStream.close();
     
   return points;
}


template<typename T, typename N>//---Prints the contents of a point set.
void Datasource2<T,N>::print(vector<N>& pointset) {
   std::cout.precision(53);
   std::cout.setf(ios::fixed);

   cout << "pointset.size=" << pointset.size() << "\n";

   for (unsigned int i=0; i<pointset.size(); i++) {
      std::cout << pointset[i].x() << " " << pointset[i].y() << "\n";
   }
}


template<typename T, typename N>//---Prints the points in Postscript.
void Datasource2<T,N>::printPS(vector<N>& pointset) {

   cout << "%!PS-Adobe-2.0 EPSF-2.0\n"
        << "1.5 1.5 scale	% scale coordinate system.\n"
        << "295 50 translate	% put origin here\n"
        << "0 0 0 setrgbcolor\n"
        << "0.3 setlinewidth\n"
        << "%----------------\n";

   // This prints out the axes with marks for every 5 units:
   cout << "newpath -500 0 moveto 500 0 lineto stroke\n"
        << "newpath 0 0 moveto 0 500 lineto stroke\n"
        << "/Helvetica findfont 3 scalefont setfont\n"
        << "-6 -8 moveto (0,0) show % Mark origo\n"
        << "0.1 setlinewidth\n";

   for (int y=5; y<305; y+=5) {
      cout << "newpath -500 " << y << " moveto 500 " << y << " lineto stroke\n"
           << "-6 " << (y-2) << " moveto (" << y << ") show\n";
   }
   for (int x=-500; x<500; x+=5) {
      cout << "newpath " << x << " 0 " << " moveto " << x << " 500 " << " lineto stroke\n";
      if (x%10 == 0) cout << x << " -6 moveto (" << x << ") show\n";
   }
   // End coordinate system.
   
   cout << "0.3 setlinewidth\n";

   // Now print out the points themselves:
   for (unsigned int i=0; i<pointset.size(); i++) {
      cout << "newpath " << pointset[i].x() << " " << pointset[i].y() << " moveto "
           << pointset[i].x() << " " << pointset[i].y() << " " << " 0.4 0 360 arc stroke\n";
   }

   cout << "%------------------------------\n"
        << "showpage	% Marks end of page";
}


/**** End of file i/o related functions. ****/

/**** Utility functions ****/

template<typename T, typename N>
void Datasource2<T,N>::translatePointset(vector<N>& pointset, float translX, float translY) {
    for(typename vector<N>::iterator it = pointset.begin(); it != pointset.end(); ++it) {
        *it = N(it->x()+translX, it->y()+translY);
    }
}

template<typename T, typename N>
void Datasource2<T,N>::dumpPointset(vector<N>& pointset) {
    for(typename vector<N>::iterator it = pointset.begin(); it != pointset.end(); ++it) {
        cout << it->x() << ", " << it->y() << "\n";
    }
    cout << "Size: " << pointset.size() << "\n";
}


/* --------------------------------------------------------------------
Re. randomPointInCircle: I found this on Mathforum.org - a description 
of how to create random points on a circle. 
See also mathworld (disk point picking).
MN oct 2008.
-----------------------------------------------------------------------
Date: 04/02/2007 at 17:41:13
From: PC
Subject: Random point within a circle.

Given a circle C with center (xc, yc) and radius R, find a random 
point within the circle.  What's the best possible way to do that in
terms of computation?



Date: 04/03/2007 at 08:57:09
From: Doctor George
Subject: Re: Random point within a circle.

Hi PC,

Thanks for writing to Doctor Math.

I assume that you are looking for the points to be uniformly
distributed in the circle.

One method is to pick random points in a square in which the circle is
circumscribed.  Then reject those points that are not also in the circle.

Another method is to select a random pair of values (M, theta) where M
is the magnitude of the vector from the center to a random point and
theta is the angle by which the vector is rotated.  The random points
then have the following form.

  x = xc + M cos(theta)
  y = yc + M sin(theta)

The distribution of theta is uniform on the interval [0, 2pi].  The
distribution of M has a ramp shape on the interval [0, R].  It is a
good exercise to prove that M has this distribution.

To compute random values of M you can use the inverse of its
cumulative distribution.  Another technique is to generate a pair of
values (r1, r2) that are uniform on the interval [0, R] and let M =
max(r1, r2).  This proof is also a good exercise.

I prefer the second method, and I typically generate values of M using
the (r1, r2) pairs.

Does that make sense?  Write again if you need more help.
*/

#endif
