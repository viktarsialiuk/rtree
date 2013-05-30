#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <vector>

#include "rtree.h"

using namespace std;
using namespace rtree;

static const int kMaxTreeElements = 1000000;

int main(int argc, char* argv[])
{
    clock_t t;
    R_tree<int> tree;
    vector<int> found;

    srand((unsigned)time(NULL));

    t = clock();
    for (int i = 0; i < kMaxTreeElements; ++i)
    {
        tree.insert(Rect(rand(), rand(), rand(), rand()), i);
    }
    t = clock() - t;
    cout << "insertion test: " << (((float)t)/CLOCKS_PER_SEC) << " seconds" << endl;



    t = clock();
    tree.rect_search(Rect(0, 0, 1000, 1000), found);
    t = clock() - t;
    cout << "rect_search test: " << (((float)t)/CLOCKS_PER_SEC) << " seconds, found " <<  found.size() << endl;

    found.clear();
    t = clock();
    tree.nearest_search(Point(0, 0), 1000, found);
    t = clock() - t;
    cout << "nearest_search test: " << (((float)t)/CLOCKS_PER_SEC) << " seconds, found " <<  found.size() << endl;

    t = clock();
    for (R_tree<int>::iterator it = tree.begin(); it != tree.end(); ++it)
    {
        *it = 0;
    }
    t = clock() - t;
    cout << "iteration test: " << (((float)t)/CLOCKS_PER_SEC) << " seconds" << endl;


    std::cin.get();
    return 0;
}