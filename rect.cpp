#include "rect.h"

#include <queue>
#include <set>
#include <algorithm>
#include <limits>
#include <math.h>

namespace rtree
{

using namespace std;



Point::Point(long x, long y)
    :
    x_(x),
    y_(y)
{
}


void
Point::set_x(long x)
{
    x_ = x;
}

long
Point::x() const
{
    return x_;
}

void
Point::set_y(long y)
{
    y_ = y;
}

long
Point::y() const
{
    return y_;
}

void
Point::shift(coord_type dx, coord_type dy)
{
    x_ += dx;
    y_ += dy;
}

void
Point::rotate(double cosx, double sinx)
{
    cosx = cosx / sqrt(cosx * cosx + sinx * sinx);
    sinx = sinx / sqrt(cosx * cosx + sinx * sinx);

    double rx = x_ * cosx - y_ * sinx;
    double ry = x_ * sinx + y_ * cosx;

    x_ = static_cast<long>(rx);
    y_ = static_cast<long>(ry);
}



Rect::Rect()
    :
    left_(numeric_limits<long>::max()),
    right_(numeric_limits<long>::min()),
    bottom_(numeric_limits<long>::max()),
    top_(numeric_limits<long>::min())
{
}

Rect::Rect(Point const& lb, Point const& rt)
    :
    left_(min(lb.x(), rt.x())),
    right_(max(lb.x(), rt.x())),
    bottom_(min(lb.y(), rt.y())),
    top_(max(lb.y(), rt.y()))
{
}

Rect::Rect(long x1, long y1, long x2, long y2)
    :
    left_(min(x1, x2)),
    right_(max(x1, x2)),
    bottom_(min(y1, y2)),
    top_(max(y1, y2))
{
}



void
Rect::set_left(long x)
{
    left_ = x;
}

long
Rect::left() const
{
    return left_;
}

void
Rect::set_right(long x)
{
    right_ = x;
}

long
Rect::right() const
{
    return right_;
}

void
Rect::set_top(long y)
{
    top_ = y;
}

long
Rect::top() const
{
    return top_;
}

void
Rect::set_bottom(long y)
{
    bottom_ = y;
}

long
Rect::bottom() const
{
    return bottom_;
}

void
Rect::set(long x1, long y1, long x2, long y2)
{
    left_ = min(x1, x2);
    right_ = max(x1, x2);
    bottom_ = min(y1, y2);
    top_ = max(y1, y2);
}


long
Rect::width() const
{
    return abs(right_ - left_);
}

long
Rect::height() const
{
    return abs(top_ - bottom_);
}


void
Rect::normalize()
{
    Rect(left_, bottom_, right_, top_).swap(*this);
}


void
Rect::rotate(double cosx, double sinx)
{
    Point lb(left_, bottom_);
    Point rt(right_, top_);

    lb.rotate(cosx, sinx);
    rt.rotate(cosx, sinx);

    Rect(lb, rt).swap(*this);
}

void
Rect::shift(long dx, long dy)
{
    left_ += dx;
    right_ += dx;

    bottom_ += dy;
    top_ += dy;
}

void
Rect::swap(Rect& oth)
{
    std::swap(left_, oth.left_);
    std::swap(bottom_, oth.bottom_);
    std::swap(right_, oth.right_);
    std::swap(top_, oth.top_);
}

double
Rect::area() const
{
    return static_cast<double>(abs(right_ - left_) * abs(top_ - bottom_));
}



double
Rect::hoverlap(Rect const& oth) const
{
    return static_cast<double>(std::min(right_, oth.right_) - std::max(left_, oth.left_));
}

double
Rect::voverlap(Rect const& oth) const
{
    return static_cast<double>(std::min(top_, oth.top_) - std::max(bottom_, oth.bottom_));
}

double
Rect::overlap(Rect const& oth) const
{
    double hover = hoverlap(oth);
    double vover = voverlap(oth);

    if (hover > 0.0 && vover > 0.0)
    {
        return hover * vover;
    }
    return 0.0;
}

void
Rect::outersect(Rect const& oth)
{
    left_ = min(left_, oth.left_);
    right_ = max(right_, oth.right_);
    bottom_ = min(bottom_, oth.bottom_);
    top_ = max(top_, oth.top_);

}


void
Rect::intersect(Rect const& oth)
{
    Rect(max(left_, oth.left_),
         max(bottom_, oth.bottom_),
         min(right_, oth.right_),
         min(top_, oth.top_)).swap(*this);
}


//void
//Rect::hintersect(Rect const& oth)
//{
//    left_ = max(left_, oth.left_);
//    right_ = min(right_, oth.right_);
//}
//
//void
//Rect::vintersect(Rect const& oth)
//{
//    bottom_ = max(bottom_, oth.bottom_);
//    top_ = min(top_, oth.top_);
//}





//return area of rect
long
area(Rect const& rect)
{
    return rect.width() * rect.height();
}


long
sqare_distance(Rect const& rect, Point const& point)
{
    long dist = 0;

    long rect_coords[] = { rect.left(), rect.right(), rect.bottom(), rect.top() };
    long point_coords[] = { point.x(), point.y() };

    for (size_t ndim = 0; ndim < sizeof(point_coords) / sizeof(long); ++ndim)
    {
        long lbound = rect_coords[2 * ndim];
        long ubound = rect_coords[2 * ndim + 1];

        if (point_coords[ndim] < lbound)
            dist += (lbound - point_coords[ndim]) * (lbound - point_coords[ndim]);
        else if (point_coords[ndim] > ubound)
            dist += (ubound - point_coords[ndim]) * (ubound - point_coords[ndim]);
    }

    return dist;
}

bool
intersection(Rect const& lhs, Rect const& rhs, Rect& res)
{
    if (!intersects(lhs, rhs))
        return false;

    res.set(max(lhs.left(), rhs.left()),
            max(lhs.bottom(), rhs.bottom()),
            min(lhs.right(), rhs.right()),
            min(lhs.top(), rhs.top()));
    return true;
}

//computes intersection of two rectangles
Rect
intersection(Rect const& lhs, Rect const& rhs)
{
    Rect tmp(lhs);
    tmp.intersect(rhs);
    return tmp;
}

void
outersection(Rect const& lhs, Rect const& rhs, Rect& res)
{
    res = lhs;
    res.outersect(rhs);
}


Rect
outersection(Rect const& lhs, Rect const& rhs)
{
    Rect tmp(lhs);
    tmp.outersect(rhs);
    return tmp;
}

//returns true if rhs intersets lhs i.e. area(intersection(lhs, rhs)) > 0
bool
intersects(Rect const& lhs, Rect const& rhs)
{
    if (lhs.left() >= rhs.right() ||
        lhs.right() <= rhs.left() ||
        lhs.bottom() >= rhs.top() ||
        lhs.top() <= rhs.bottom())
        return false;

    return true;
}

//returns true if circle with origin point and radius intersects lhs
bool
intersects(Rect const& rect, Point const& point, long radius)
{
    return sqare_distance(rect, point) < radius * radius;
}





class UnionFind
{
public:
    UnionFind(size_t size)
        :
        parent_(size, 0),
        weight_(size, 1)
    {
        for (size_t i = 0; i < parent_.size(); ++i)
            parent_[i] = i;
    }

    bool connected(size_t i, size_t j) const
    {
        return root(i) == root(j);
    }

    //union
    void connect(size_t i, size_t j)
    {
        i = root(i);
        j = root(j);

        if (i == j)
            return;

        if (weight_[i] >= weight_[j])
        {
            parent_[j] = i;
            weight_[i] += weight_[j];
        }
        else
        {
            parent_[i] = j;
            weight_[j] += weight_[i];
        }
    }

    size_t root(size_t i) const
    {
        while (i != parent_[i])
            i = parent_[i];
        return i;
    }
private:
    vector<size_t> parent_;
    vector<size_t> weight_;
};




struct QueueEvent
{
    QueueEvent(size_t ind, long coord, bool insert)
        :
        ind_(ind),
        coord_(coord),
        insert_(insert)
    {
    }

    bool operator<(QueueEvent const& oth) const
    {
        return coord_ < oth.coord_;
    }

    size_t      ind_;
    long        coord_;
    bool        insert_;
};

void merge(Rect_vector const& src, Rect_vector& dest)
{
    deque<QueueEvent>   queue;
    Rect_vector         res(src.size());
    UnionFind           ufind(src.size());
    set<size_t>         active;

    //initialization
    for (size_t i = 0; i < src.size(); ++i)
    {
        queue.push_back(QueueEvent(i, src[i].left(), true));
        queue.push_back(QueueEvent(i, src[i].right(), false));
    }

    //stable to be sure that insert event is before remove
    stable_sort(queue.begin(), queue.end());
    while (!queue.empty())
    {
        QueueEvent cur = queue.front();
        queue.pop_front();

        if (cur.insert_)
        {
            //process insert event
            //connect overlaped rectangles

            set<size_t>::iterator it;
            for (it = active.begin(); it != active.end(); ++it)
            {
                //we found active rectangle in queue that overlaps current rectangle
                //remove it from search
                if (src[*it].overlap(src[cur.ind_]))
                    ufind.connect(*it, cur.ind_);
            }

            active.insert(cur.ind_);
        }
        else
        {
            //process remove event
            active.erase(cur.ind_);
        }
    }


    for (size_t i = 0; i < res.size(); ++i)
        res[ufind.root(i)].outersect(src[i]);

    //copy to destination
    dest.clear();
    for (size_t i = 0; i < res.size(); ++i)
        if (ufind.root(i) == i)
            dest.push_back(res[i]);
}




}//namespace rtree