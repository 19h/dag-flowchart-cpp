#ifndef POINTF_H
#define POINTF_H

#include <vector>

/**
 * @brief Simple 2D point with double precision coordinates.
 * 
 * A lightweight replacement for QPointF that provides basic point operations
 * needed for graph layout calculations.
 */
struct PointF {
    double x;
    double y;

    PointF() : x(0), y(0) {}
    PointF(double x, double y) : x(x), y(y) {}

    void setX(double newX) { x = newX; }
    void setY(double newY) { y = newY; }

    // Allow in-place modification like Qt's rx()/ry()
    double& rx() { return x; }
    double& ry() { return y; }
    const double& rx() const { return x; }
    const double& ry() const { return y; }

    PointF& operator-=(const PointF& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }
};

/// A polyline is a sequence of points forming a path
using Polyline = std::vector<PointF>;

#endif // POINTF_H
