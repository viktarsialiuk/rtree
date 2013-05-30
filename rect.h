#pragma once
#include <vector>
#include <stack>
#include <functional>

namespace azone
{

class Point
{
public:
    typedef long            coord_type;


    Point(long x = 0L, long y = 0L);

    void        set_x(long x);
    long        x() const;

    void        set_y(long y);
    long        y() const;

    void        shift  (coord_type dx, coord_type dy);
    void        rotate (double cosx, double sinx);

private:
    long    x_;
    long    y_;
};


struct Point_compare_x
{
    int compare(Point const& lhs, Point const& rhs) const
    {
        if (lhs.x() == rhs.x())
        {
            return 0;
        }
        return lhs.x() < rhs.x() ? -1 : 1;
    }

    bool operator()(Point const& lhs, Point const& rhs) const
    {
        return compare(lhs, rhs) < 0;
    }
};



struct Point_compare_y
{
    int compare(Point const& lhs, Point const& rhs) const
    {
        if (lhs.y() == rhs.y())
        {
            return 0;
        }
        return lhs.y() < rhs.y() ? -1 : 1;
    }

    bool operator()(Point const& lhs, Point const& rhs) const
    {
        return compare(lhs, rhs) < 0;
    }
};



class Rect
{
public:
    typedef Point::coord_type          coord_type;

    Rect();
    Rect(Point const& lb, Point const& rt);
    Rect(long x1, long y1, long x2, long y2);
    Rect(double x1, double y1, double x2, double y2);

    void        set_left    (long x);
    long        left        () const;

    void        set_right    (long x);
    long        right        () const;

    void        set_top        (long y);
    long        top            () const;

    void        set_bottom    (long y);
    long        bottom        () const;

    void        set            (long x1 = 0L, long y1 = 0L, long x2 = 0L, long y2 = 0L);
    void        set            (double x1 = 0.0, double y1 = 0.0, double x2 = 0.0, double y2 = 0.0);
    //Point        left_bottom() const;
    //Point        right_top() const;

    Point       left_bottom () const
    {
        Point t(left_, bottom_);
        return t;
    }

    Point       right_top () const
    {
        Point t(right_, top_);
        return t;
    }


    void        rotate      (double cosx, double sinx);
    void        shift       (long dx, long dy);

    long        width        () const;
    long        height        () const;

    double        area        () const;
    double        hoverlap    (Rect const& oth) const;
    double        voverlap    (Rect const& oth) const;
    double        overlap        (Rect const& oth) const;

    //merge two rectangle into one bounding box
    void        outersect    (Rect const& oth);
    //intersection of two rectangles
    void        intersect    (Rect const& oth);
    //void        hintersect   (Rect const& oth);
    //void        vintersect   (Rect const& oth);

    void        normalize();
    void        swap(Rect& oth);

private:
    long    left_;
    long    bottom_;
    long    right_;
    long    top_;
};

typedef std::vector<Rect>           Rect_vector;
typedef std::stack<Rect>            Rect_stack;


//returns area of rect
long    area            (Rect const& rect);

//returns distance from point to rect
long    sqare_distance  (Rect const& rect, Point const& point);

//returns true if rhs intersets lhs i.e. area(intersection(lhs, rhs)) > 0
bool    intersects      (Rect const& lhs, Rect const& rhs);

//returns true if circle with origin point and radius intersects lhs
bool    intersects      (Rect const& rect, Point const& point, long radius);

//computes intersection of two rectangles
bool    intersection    (Rect const& lhs, Rect const& rhs, Rect& res);

//computes intersection of two rectangles
Rect    intersection    (Rect const& lhs, Rect const& rhs);


//computes bounding box of two rectangles
void    outersection    (Rect const& lhs, Rect const& rhs, Rect& res);
Rect    outersection    (Rect const& lhs, Rect const& rhs);


struct Rect_compare_left : public std::binary_function<Rect, Rect, int>
{
    int compare(Rect const& lhs, Rect const& rhs) const
    {
        //return lhs.left() - rhs.left()
        if (lhs.left() == rhs.left())
        {
            return 0;
        }
        return lhs.left() < rhs.left() ? -1 : 1;
    }

    int operator()(Rect const& lhs, Rect const& rhs) const
    {
        return compare(lhs, rhs);
    }
};

struct Rect_compare_right : public std::binary_function<Rect, Rect, int>
{
    int compare(Rect const& lhs, Rect const& rhs) const
    {
        if (lhs.right() == rhs.right())
        {
            return 0;
        }
        return lhs.right() < rhs.right() ? -1 : 1;
    }

    int operator()(Rect const& lhs, Rect const& rhs) const
    {
        return compare(lhs, rhs);
    }
};

struct Rect_compare_top : public std::binary_function<Rect, Rect, int>
{
    int compare(Rect const& lhs, Rect const& rhs) const
    {
        if (lhs.top() == rhs.top())
        {
            return 0;
        }
        return lhs.top() < rhs.top() ? -1 : 1;
    }

    int operator()(Rect const& lhs, Rect const& rhs) const
    {
        return compare(lhs, rhs);
    }
};


struct Rect_compare_bottom : public std::binary_function<Rect, Rect, int>
{
    int compare(Rect const& lhs, Rect const& rhs) const
    {
        if (lhs.bottom() == rhs.bottom())
        {
            return 0;
        }
        return lhs.bottom() < rhs.bottom() ? -1 : 1;
    }

    int operator()(Rect const& lhs, Rect const& rhs) const
    {
        return compare(lhs, rhs);
    }
};

template<class T>
struct Comparator_asc
{
    Comparator_asc(T comparator) : comparator_(comparator) {}
    bool operator()(typename T::first_argument_type lhs, typename T::second_argument_type rhs) const
    {
        return comparator_(lhs, rhs) < 0;
    }

private:
    T comparator_;
};

template<class T>
struct Comparator_desc
{
    Comparator_desc(T comparator) : comparator_(comparator) {}
    bool operator()(typename T::first_argument_type lhs, typename T::second_argument_type rhs) const
    {
        return comparator_(lhs, rhs) > 0;
    }

private:
    T comparator_;
};


template<class T>
struct Comparator_equal
{
    Comparator_equal(T comparator) : comparator_(comparator) {}
    bool operator()(typename T::first_argument_type lhs, typename T::second_argument_type rhs) const
    {
        return comparator_(lhs, rhs) == 0;
    }

private:
    T comparator_;
};

template<class T>
Comparator_asc<T> compare_asc(T comparator)
{
    return Comparator_asc<T>(comparator);
}

template<class T>
Comparator_desc<T> compare_desc(T comparator)
{
    return Comparator_desc<T>(comparator);
}

template<class T>
Comparator_equal<T> compare_equal(T comparator)
{
    return Comparator_equal<T>(comparator);
}


void merge(Rect_vector const& src, Rect_vector& dest);


}//namespace azone