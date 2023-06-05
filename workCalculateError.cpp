#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using namespace std;

// This function returns the distnace **in meters** between to locations
// specified in terms of lat and lon.
double GPSDistance(double latC1, double lonC1, double latR2, double lonR2) {
    double m_per_deg_lat = 111079;
    double m_per_deg_lon = 82463;
    double latDifference = (latC1 - latR2) * m_per_deg_lat;
    double lonDifference = (lonC1 - lonR2) * m_per_deg_lon;

    double distance = sqrt(pow(latDifference, 2) + pow(lonDifference, 2));

    return distance;
}

double GPSAngle(double latC1, double lonC1, double latR2, double lonR2) {

    double changeX = lonR2 - lonC1;
    double changeY = latR2 - latC1;
    double result = -atan2(changeY, changeX);

    return result;
}

int closestRefrenceIndex(double latC1, double lonC1, double arrRLat[],
                         double arrRLon[], long n) {
    double closestDistance = DBL_MAX;
    int closestIndex = -1;

    for (int i = 0; i < n; i++) {
        double RLat = arrRLat[i];
        double Rlon = arrRLon[i];
        double tempClosest = GPSDistance(latC1, lonC1, RLat, Rlon);
        if (tempClosest < closestDistance) {
            double latR2point = arrRLat[i];
            double lonR2point = arrRLon[i];
            closestDistance = tempClosest;
            closestIndex = i;
        }
    }
    return closestIndex;
}

typedef pair<double, double> Point;

Point GetClosestPoint(Point A, Point B, Point P) {
    Point A_to_P;
    Point A_to_B;
    Point closestPoint;

    // thinkg of ".first" as like ".x" and ".second" as like ".y"
    A_to_P.first = P.first - A.first;
    A_to_P.second = P.second - A.second;

    A_to_B.first = B.first - B.first;
    A_to_B.second = B.second - A.second;

    double squaredMagnitudeof_A_to_B =
        pow(A_to_B.first, 2) + pow(A_to_B.second, 2);

    double dotProductof_AtoB_AtoP =
        (A_to_P.first * A_to_B.first) + (A_to_P.second * A_to_B.second);

    double normalizedDistanceFrom_A_to_closestPoint =
        dotProductof_AtoB_AtoP / squaredMagnitudeof_A_to_B;

    closestPoint.first =
        A.first + A_to_B.first * normalizedDistanceFrom_A_to_closestPoint;
    closestPoint.second =
        A.second + A_to_B.second * normalizedDistanceFrom_A_to_closestPoint;

    return closestPoint;
}

// // This function computes the distance **in meters** between the two points
// // specified in terms of latitude and longitude.
// double findDistanceError(Point P, Point closestPoint) {
//     // Use GPSDistance
//     double errorDistance = pow(P.first - closestPoint.first, 2) +
//                            pow(P.second - closestPoint.second, 2);

//     return pow(errorDistance, 0.5);
// }

int main() {

    double arrCLat[10] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    double arrCLon[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double arrRLat[10] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1};
    double arrRLon[10] = {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1, 10.1};

    double errors[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int n = 10;

    for (int c = 0; c < n; c++) {
        int secondClosestIndex = -1;
        // Find index of closest R point
        int indexClosestR =
            closestRefrenceIndex(arrCLat[c], arrCLon[c], arrRLat, arrCLat, n);

        // Look on either side (before and after) that point to find NEXT
        // closest point

        if (indexClosestR == 0) {
            secondClosestIndex = indexClosestR + 1;
        } else if (indexClosestR == n - 1) {
            secondClosestIndex = indexClosestR - 1;

        } else {
            double distanceBefore =
                GPSDistance(arrCLat[c], arrCLon[c], arrRLat[indexClosestR - 1],
                            arrCLon[indexClosestR - 1]);
            double distanceAfter =
                GPSDistance(arrCLat[c], arrCLon[c], arrRLat[indexClosestR + 1],
                            arrCLon[indexClosestR + 1]);

            if (distanceBefore < distanceAfter) {
                secondClosestIndex = indexClosestR - 1;
            } else {
                secondClosestIndex = indexClosestR + 1;
            }
        }

        // Find angle between two reference points. Find angle
        //  between comparison point and reference line (+/- 90 degrees of
        //  previous angle)
        // GPSAngle(arrCLat[c], arrCLon[c], arrRLat[indexClosestR],
        //          arrRLon[indexClosestR]);
        // Find point on reference line closest to comparison
        //  point
        Point A, B, P;
        A.first = arrRLat[indexClosestR];
        A.second = arrRLon[indexClosestR];
        B.first = arrRLat[secondClosestIndex];
        B.second = arrRLon[secondClosestIndex];
        P.first = arrCLat[c];
        P.second = arrCLon[c];
        Point closestPoint = GetClosestPoint(A, B, P);
        // Find distance between last point and comparison point...this
        //  is the error. Set error for this point.
        // comparison point(c1) - refrence point(r1);
        double distanceError = GPSDistance(
            P.first, P.second, closestPoint.first, closestPoint.second);
        errors[c] = distanceError;
        // cout << distanceError;
    }

    return 0;
}
