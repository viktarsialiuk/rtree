#pragma once
#include "rect.h"

#include <math.h>           //abs
#include <assert.h>         //assert
#include <algorithm>        //swap
#include <functional>       //std::equal_to
#include <iterator>         //std::iterator
#include <memory>           //std::allocator
#include <limits>           //std::numeric_limits

//#include <boost/aligned_storage.hpp>
//#include <boost/type_traits.hpp>
#include <type_traits>

namespace rtree
{

namespace details
{

template<class T>
struct lowest
{
    static T value()
    {
        return std::numeric_limits<T>::min();
    }
};

template<>
struct lowest<float>
{
    static float value()
    {
        return -std::numeric_limits<float>::max();
    }
};

template<>
struct lowest<double>
{
    static double value()
    {
        return -std::numeric_limits<double>::max();
    }
};


template<>
struct lowest<long double>
{
    static long double value()
    {
        return -std::numeric_limits<long double>::max();
    }
};



template<class T, size_t F>
struct tree_node
{
    typedef Rect                                                    rect_type;
    typedef T                                                       value_type;

    tree_node(size_t level = 0, rect_type const& bbox = rect_type())
        : level_(level), bbox_(bbox), count_(0), parent_(0), index_(0)
    {
#ifdef RTREE_DEBUG
        ++allocated;
#endif
    }


    ~tree_node()
    {
#ifdef RTREE_DEBUG
        --allocated;
#endif
    }

    bool is_leaf() const
    {
        return 0 == level_;
    }

    bool empty() const
    {
        return 0 == count_;
    }

    size_t level() const
    {
        return level_;
    }
    size_t count() const
    {
        return count_;
    }

    size_t index() const
    {
        return index_;
    }
    tree_node* parent() const
    {
        return parent_;
    }

    rect_type& bbox()
    {
        return bbox_;
    }

    rect_type const& bbox() const
    {
        return bbox_;
    }

    //NOTE: data is constructed only for leaf nodes
    //that's why it can not be accessed in high level nodes
    value_type* data()
    {
        assert(level_ == 0);
        return static_cast<value_type*>(static_cast<void*>(&data_));
    }

    void clear()
    {
        //clear
        count_ = 0;
    }

    void reset_links(tree_node* parent, size_t index)
    {
        parent_ = parent;
        index_ = index;
    }

    //erase child node
    tree_node* erase_child(size_t ind)
    {
        assert(ind < count_);

        tree_node* res = childs_[ind];
        res->reset_links(NULL, 0);//reset parent child links

        //shift all child nodes
        for (size_t i = ind; i + 1 < count_; ++i)
        {
            childs_[i] = childs_[i + 1];
            childs_[i]->index_ = i;
        }

        --count_;
        return res;
    }

    //returns child branch
    tree_node* get_child(size_t ind)
    {
        assert(ind < count_);
        return childs_[ind];
    }

    tree_node const* get_child(size_t ind) const
    {
        assert(ind < count_);
        return childs_[ind];
    }

    //insert new child branch into the node and reset child/parent links
    bool add_child(tree_node* child)
    {
        if (count_ >= F)
        {
            return false;
        }

        //reset parent child links
        child->reset_links(this, count_);
        childs_[count_++] = child;
        return true;
    }


private:
    size_t                              level_;         //level of the node
    size_t                              count_;         //first unused branch
    size_t                              index_;         //index of current branch in parent node (required for iterator implementation)
    tree_node*                          parent_;        //parent node (required for iterator implementation)
    rect_type                           bbox_;          //bounding box
    tree_node*                          childs_[F];

    typename std::aligned_storage<sizeof(T),
                                  std::alignment_of<T>::value >::type data_;

#ifdef RTREE_DEBUG
public:
    static size_t allocated;
#endif


private:
    //disallow copy and assign
    tree_node(tree_node const&);
    tree_node& operator=(tree_node const&);
};

#ifdef RTREE_DEBUG
template<class T, size_t F>
size_t tree_node<T, F>::allocated = 0;
#endif



//tree iterator implementation
template<class T, size_t F>
struct tree_iterator
{
    typedef Rect                                                    rect_type;
    typedef tree_node<typename std::remove_const<T>::type, F>       tree_node;
    typedef T                                                       value_type;

    //iterator required definitions
    typedef std::forward_iterator_tag                               iterator_category;
    typedef ptrdiff_t                                               difference_type;
    typedef ptrdiff_t                                               distance_type;
    typedef value_type*                                             pointer;
    typedef value_type&                                             reference;

public:
    tree_iterator(tree_node* root = NULL)
        : node_(move_to_child(root))
    {
    }

    tree_iterator(tree_iterator< typename std::remove_const<T>::type, F > const& oth)//allow constructing from non-const iterator
        :
        node_(oth.node())
    {
    }

    value_type& operator*() const
    {
        assert(node_ && node_->level() == 0);
        return *(node_->data());
    }

    value_type* operator->() const
    {
        assert(node_ && node_->level() == 0);
        return node_->data();
    }

    tree_iterator& operator++()
    {
        increment();
        return *this;
    }

    tree_iterator operator++(int)
    {
        tree_iterator temp(*this);
        ++*this;
        return temp;
    }

    void swap(tree_iterator& oth)
    {
        std::swap(node_, oth.node_);
    }

    bool operator==(tree_iterator const& oth) const
    {
        return node_ == oth.node_;
    }
    bool operator!=(tree_iterator const& oth) const
    {
        return !this->operator==(oth);
    }


    tree_node* node() const
    {
        return node_;
    }

private:
    void increment()
    {
        if (!node_)
        {
            //end
            return;
        }

        node_ = move_to_parent(node_);
        node_ = move_to_child(node_);
    }

    //move current position to the first left leaf of the node
    tree_node* move_to_parent(tree_node* leaf)
    {
        //move to parent and skip all visited nodes
        for (tree_node* node = leaf; node && node->parent(); node = node->parent())
        {
            if (node->index() + 1 < node->parent()->count())//check if all nodes in parent node are already visited
            {
                //return next unvisited child
                return node->parent()->get_child(node->index() + 1);
            }
        }

        //all nodes are already visited
        return NULL;
    }

    //move current position to the first left leaf of the node
    tree_node* move_to_child(tree_node* node)
    {
        //move to first left leaf
        while (node && node->level() > 0)
        {
            node = node->empty() ? NULL : node->get_child(0);
        }

        return node;
    }

private:
    tree_node* node_;
};





//struct for storing data found during nearest neighbour search
//and sorting branches by distance
template<class T>
struct heap_data
{
    heap_data(Point::coord_type distance = 0, T const& data = T())
        : distance_(distance), data_(data) {}

    bool operator<(heap_data const & oth)
    {
        return distance_ < oth.distance_;
    }

    T data() const
    {
        return data_;
    }

    Point::coord_type distance() const
    {
        return distance_;
    }

private:
    T                   data_;
    Point::coord_type   distance_;
};


template<class T>
class K_nearest_heap
{
public:
    K_nearest_heap(size_t max_size) : max_size_(max_size) {}

    bool full() const
    {
        return max_size_ <= data_.size();
    }

    bool empty() const
    {
        return data_.empty();
    }

    Point::coord_type max_distance() const
    {
        if (full())
        {
            return data_.front().distance();
        }

        return std::numeric_limits<Point::coord_type>::max();
    }

    void pop()
    {
        std::pop_heap(data_.begin(), data_.end());
        data_.pop_back();
    }


    bool push(heap_data<T> const& item)
    {
        if (max_distance() < item.distance())
        {
            //heap is full and current max distance is lower than item distance - skip current item
            return false;
        }

        if (full())
        {
            pop();
        }

        data_.push_back(item);
        std::push_heap(data_.begin(), data_.end());
        return true;
    }

    heap_data<T> const& top()
    {
        return data_.front();
    }

private:
    size_t                          max_size_;
    std::vector< heap_data<T> >     data_;
};


}//namespace details


//R_tree implementation
template<class T, size_t F = 2, class Alloc = std::allocator<T> >
class R_tree
{

private:
    enum
    {
        kMaxFill = 2 * F + 1,
        kMinFill = F,
        kMinRootFill = 2
    };


public:
    typedef T                                                       value_type;
    typedef Rect                                                    rect_type;
    typedef Point                                                   point_type;
    typedef Point::coord_type                                       coord_type;

    //typedef std::pair<rect_type, data_type>                         value_type;
    typedef value_type&                                             reference;
    typedef value_type*                                             pointer;
    typedef value_type const&                                       const_reference;
    typedef value_type const*                                       const_pointer;

    typedef details::tree_node<value_type, kMaxFill>                tree_node;

    typedef details::tree_iterator<value_type, kMaxFill>            iterator;
    typedef details::tree_iterator<const value_type, kMaxFill>      const_iterator;

    //R_tree definition
    R_tree()
        :
        cached_size_(0),
        root_(NULL)
    {
    }

    ~R_tree()
    {
        clean(root_);
    }

    iterator begin()
    {
        return iterator(root_);
    }
    const_iterator begin() const
    {
        return const_iterator(root_);
    }

    iterator end()
    {
        return iterator(NULL);
    }
    const_iterator end() const
    {
        return const_iterator(NULL);
    }

    iterator find(rect_type const& query, value_type const& val) const
    {
        if (!root_)
            return iterator(NULL);

        return find_internal(root_, query, std::bind1st(std::equal_to<value_type>(), val));
    }

    template<class Pred>
    iterator find(rect_type const& query, Pred const& pred) const
    {
        if (!root_)
            return iterator(NULL);

        return find_internal(root_, query, pred);
    }

    void erase(const_iterator wh)
    {
        if (wh != end())
        {
            erase_internal(wh.node());
            --cached_size_;
        }
    }

    void insert(rect_type const& bbox, value_type const& data)
    {
        //create new leaf node
        tree_node* leaf = create_node(bbox, data);

        try
        {
            insert_internal(leaf);
        }
        catch(bad_alloc&)
        {
            //cached_size_ = compute_size();
            throw;
        }

        //update cached size
        ++cached_size_;
    }


    void clear()
    {
        R_tree().swap(*this);
    }

    size_t size() const
    {
        return cached_size_;
    }

    size_t compute_size() const
    {
        return std::distance(begin(), end());
    }

    bool empty() const
    {
        return 0 == cached_size_;
    }

    template<class Container>
    void rect_search2(rect_type const& query, Container& found) const
    {
        if (!root_)
            return;

        rect_search2(root_, query, found);
    }

    //returns all elements that intersect query rectnagle
    template<class Container>
    void rect_search(rect_type const& query, Container& found) const
    {
        if (!root_)
            return;

        rect_search(root_, query, found);
    }

    //returns all elemens that intersect circle with origin
    template<class Container>
    void rad_search(point_type const& origin, coord_type dist, Container& found) const
    {
        if (!root_)
            return;

        rad_search(root_, origin, dist, found);
    }

    template<class Container>
    void nearest_search(point_type const& query, size_t k, Container& found) const
    {
        if (!root_)
            return;

        details::K_nearest_heap<value_type> heap(k);
        //nearest search using heap
        nearest_search(root_, query, heap);

        //copy to output container
        while (!heap.empty())
        {
            found.push_back(heap.top().data());
            heap.pop();
        }
    }

    void swap(R_tree& oth)
    {
        std::swap(cached_size_, oth.cached_size_);
        std::swap(root_, oth.root_);
        std::swap(nalloc_, oth.nalloc_);
    }

private:

    void destroy_node(tree_node* node)
    {
        if (node->is_leaf())
            node->data()->~value_type();    //destroy value for leaf nodes

        node->~tree_node();
        nalloc_.deallocate(node, 1);
    }

    //used to create non-leaf nodes
    tree_node* create_node(size_t level)
    {
        assert(level > 0);

        tree_node* node = nalloc_.allocate(1);
        new (static_cast<void*>(node)) tree_node(level);//no throw
        return node;
    }

    //used to create leaf nodes
    tree_node* create_node(rect_type const& bbox, value_type const& val)
    {
        tree_node* node = nalloc_.allocate(1);
        new (static_cast<void*>(node)) tree_node(0, bbox);//no throw

        try
        {
            new (static_cast<void*>(node->data())) value_type(val);
        }
        catch(...)
        {
            node->~tree_node();
            nalloc_.deallocate(node, 1);
            throw;
        }

        return node;
    }

    //predicates
    static bool clean_non_leaf(tree_node* node)
    {
        return !node->is_leaf();
    }
    static bool clean_all(tree_node* node)
    {
        return true;
    }
    void clean(tree_node* node, bool (*pred)(tree_node*) = clean_all)
    {
        if (!node || !pred(node))
        {
            return;
        }

        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            clean(node->get_child(ind), pred);
        }

        destroy_node(node);
    }

    template<class Pred>
    iterator find_internal(tree_node* node, rect_type const& query, Pred const& pred) const
    {
        assert(node);

        if (!intersects(node->bbox(), query))
            return iterator(NULL);
        else if (node->is_leaf())
        {
            if (pred(*(node->data())))
                return iterator(node);

            return iterator(NULL);
        }

        //iterate untill nothing is found
        iterator found, end;
        for (size_t ind = 0; ind < node->count() && found == end; ++ind)
        {
            found = find_internal(node->get_child(ind), query, pred);
        }
        return found;
    }

    template<class Container>
    void rect_search2(tree_node* node, rect_type const& query, Container& found) const
    {
        assert(node);

        if (!intersects(node->bbox(), query))
            return;
        else if (node->is_leaf())
        {
            found.push_back(iterator(node));
            return;
        }

        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            //search recursively
            rect_search2(node->get_child(ind), query, found);
        }
    }

    template<class Container>
    void rect_search(tree_node* node, rect_type const& query, Container& found) const
    {
        assert(node);

        if (!intersects(node->bbox(), query))
            return;
        else if (node->is_leaf())
        {
            found.push_back(*(node->data()));
            return;
        }

        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            rect_search(node->get_child(ind), query, found);
        }
    }

    template<class Container>
    void rad_search(tree_node* node, point_type const& origin, coord_type dist, Container& found) const
    {
        assert(node);

        if (!intersects(node->bbox(), origin, dist))
            return;
        else if (node->is_leaf())
        {
            found.push_back(*(node->data()));
            return;
        }

        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            rad_search(node->get_child(ind), origin, dist, found);
        }
    }

    void nearest_search(tree_node* node, point_type const& query, details::K_nearest_heap<value_type>& heap) const
    {
        assert(node);

        if (node->is_leaf())
        {
            coord_type current = sqare_distance(node->bbox(), query);
            heap.push(details::heap_data<value_type>(current, *(node->data())));

            return;
        }

        //sort bracnhes by distance to query point
        std::vector< details::heap_data<size_t> > sorted;
        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            coord_type distance = sqare_distance(node->get_child(ind)->bbox(), query);
            sorted.push_back(details::heap_data<size_t>(distance, ind));
        }
        std::sort(sorted.begin(), sorted.end());


        for (size_t ind = 0; ind < sorted.size(); ++ind)
        {
            size_t      current_branch = sorted[ind].data();
            coord_type  current_distance = sorted[ind].distance();

            //check if heap contains k elemenents
            if (heap.max_distance() < current_distance)
            {
                //do not search anymore since current distance is bigger than the worst distance in the heap
                break;
            }

            //otherwise: search child nodes recursively
            nearest_search(node->get_child(current_branch), query, heap);
        }
    }

    void insert_internal(tree_node* leaf)
    {
        assert(leaf && leaf->level() == 0);

        tree_node* pnode = NULL;
        if (root_ && root_->level() > 0)
        {
            pnode = root_;
            while (pnode->level() > 1)
            {
                size_t ind = pick_branch(pnode, leaf->bbox());
                pnode = pnode->get_child(ind);
            }
        }

        tree_node* split = leaf;
        while (pnode)
        {
            //update parent bbox
            //if pnode is being splitted then we recompute bbox in quadratic split
            //if it isn't then bbox has to be outersected with child bbox
            outersection(pnode->bbox(), leaf->bbox(), pnode->bbox());
            if (split)
            {
                assert(validate_bbox(split));

                split = add_branch(pnode, split); //add split node to parent and update parent bbox
            }

            assert(validate_bbox(pnode));

            pnode = pnode->parent();
        }

        if (split)
        {
            //root_ is NULL or was splitted
            if (root_)
                root_ = split_root(split);
            else
            {
                assert(split == leaf);

                root_ = split;
            }

            assert(validate_bbox(root_));
        }
    }

    void erase_internal(tree_node* leaf)
    {
        //we can't delete leaf from empty tree
        assert(root_ && leaf && leaf->level() == 0);

        //move to first parent with childs count > KMinFill
        size_t      index = leaf->index();
        tree_node*  pnode = leaf->parent();

        while (pnode && !is_filled_node(pnode))
        {
            index = pnode->index();
            pnode = pnode->parent();
        }

        tree_node* erased;
        if (!pnode)
        {
            //erase root_ node
            erased = root_;
            root_ = NULL;
        }
        else
        {
            //erase child node in pnode
            erased = pnode->erase_child(index);
        }

        //recompute bboxes for all parents nodes
        while (pnode)
        {
            pnode->bbox() = recompute_bbox(pnode);
            pnode = pnode->parent();
        }

        //re-insert all child leaf nodes
        iterator it(erased);
        while(it != end())
        {
            tree_node* cur = (it++).node();
            if (cur == leaf)
                continue;

            //note that cur node is being re-inserted
            //that's why cur->parent_ and cur->index are updated but it doesn't affect other nodes in erased sub tree
            //that's why iterator is still valid and it points to the next leaf node
            cur->reset_links(NULL, 0);
            insert_internal(cur);
       }

        //we do not need to destroy leaf nodes because we have re-inserted them
        clean(erased, clean_non_leaf);
        destroy_node(leaf);
        return;
    }

    tree_node* split_root(tree_node* split)
    {
        assert(root_ && split && split->level() == root_->level());

        tree_node* nroot;
        try
        {
            nroot = create_node(root_->level() + 1);
        }
        catch (std::bad_alloc&)
        {
            clean(split, clean_all);
            throw;
        }

        nroot->add_child(root_);
        nroot->add_child(split);
        nroot->bbox() = recompute_bbox(nroot);
        return nroot;
    }

    rect_type recompute_bbox(tree_node const* node)
    {
        assert(node && node->level() > 0);

        rect_type bbox;

        for (size_t i = 0; i < node->count(); ++i)
            outersection(bbox, node->get_child(i)->bbox(), bbox);

        return bbox;
    }

    //returns true if node fill is > than MinFill
    bool is_filled_node(tree_node* node)
    {
        assert(node);

        if (node->level() == 0)
            return false;
        if (node == root_)
            return node->count() > kMinRootFill;

        return node->count() > kMinFill;
    }
    


    /*
    [Pick branch]
        Select branch which needs least enlargement to include new branch. When there are more
        qualify entries, the entry with the rectangle of the smallest area is chosen
    */

    size_t pick_branch(tree_node* node, rect_type const& bbox) const
    {
        assert(node && node->count() > 0);

        size_t      res = 0;
        coord_type  min_area = 0;
        coord_type  increase = std::numeric_limits<coord_type>::max();
        rect_type   merged;

        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            //compute bounding box
            outersection(bbox, node->get_child(ind)->bbox(), merged);
            //current branch area
            coord_type branch_area = area(node->get_child(ind)->bbox());

            //required enlargement
            coord_type increase_area = area(merged) - branch_area;

            //choose branch with smallest enlargement and smallest area
            if (increase_area < increase)
            {
                res = ind;
                increase = increase_area;
                min_area = branch_area;
            }
            else if (increase_area == increase &&
                     branch_area < min_area)
            {
                res = ind;
                min_area = branch_area;
            }
        }

        return res;
    }

    /*
    Quadratic Split pseudo-code:
        Divide a set of kMaxFill + 1 index entries into two groups.
        QS1 [Pick first entry for each group] start PickSeeds to find two entries to be
            the first elements of the groups. Assign each to a group.
        QS2 [Check if done] If all entries have been assigned, then break up. If a group
            has too few entries that all the rest must be assigned to it in order for it to
            have the minimum number m, then assign them and break up
    */
    void quadractic_split(tree_node* node, tree_node* split, tree_node* child) /*throw()*/
    {
        assert(child && node && node->count() == kMaxFill);

        tree_node* branches[kMaxFill + 1];
        size_t groups[kMaxFill + 1];

        const size_t branches_size = node->count() + 1;

        branches[0] = child;
        for (size_t ind = 0; ind < node->count(); ++ind)
        {
            branches[ind + 1] = node->get_child(ind);
        }

        size_t group_0 = 0;
        size_t group_1 = 0;

        //find groups
        pick_seeds(branches, branches_size, group_0, group_1);

        //split branches on two groups
        distribute_branches(branches, branches_size, groups, group_0, group_1);

        //assign nodes
        node->clear();
        for (size_t ind = 0; ind < branches_size; ++ind)
        {
            if (groups[ind] == group_0 + 1)
            {
                node->add_child(branches[ind]);
            }
            else
            {
                split->add_child(branches[ind]);
            }
        }
        node->bbox() = recompute_bbox(node);
        split->bbox() = recompute_bbox(split);
    }

    /*
    QS3 [Select entry to assign] start PickNext to choose the next entry to assign.
        This is put in the group whose area has to be least enlarged. If the algorithm
        enlarges both groups with the same size, then add to the group with smaller
        area, then to the one with fewer entries, then to either. Repeat from QS2.
    */
    void distribute_branches(tree_node** branches, size_t count,
                             size_t* groups,
                             size_t group_0, size_t group_1)
    {
        //clear groups
        memset(groups, 0, sizeof(size_t) * count);

        groups[group_0] = group_0 + 1;
        groups[group_1] = group_1 + 1;

        //number element in groups
        size_t      nelems[] = { 1, 1 };
        size_t      groupids[] = { group_0, group_1 };
        rect_type   bboxes[] = { branches[group_0]->bbox(), branches[group_1]->bbox() };

        rect_type merged;

        while (nelems[0] + nelems[1] < count &&
               nelems[0] < count / 2 &&
               nelems[1] < count / 2)
        {
            size_t          best_branch = pick_next(branches, count, groups, bboxes[0], bboxes[1]);
            coord_type      increase_0 = 0;
            coord_type      increase_1 = 0;
            coord_type      area_0 = area(bboxes[0]);
            coord_type      area_1 = area(bboxes[1]);

            outersection(branches[best_branch]->bbox(), bboxes[0], merged);
            increase_0 = area(merged) - area_0;

            outersection(branches[best_branch]->bbox(), bboxes[1], merged);
            increase_1 = area(merged) - area_1;


            size_t group = 0;
            if (increase_0 != increase_1)
            {
                group = increase_0 < increase_1 ? 0 : 1;
            }
            else if (area_0 != area_1)
            {
                //the same increase
                group = area_0 < area_1 ? 0 : 1;
            }
            else
            {
                group = nelems[0] < nelems[1] ? 0 : 1;
            }

            ++nelems[group];
            groups[best_branch] = groupids[group] + 1;
            outersection(bboxes[group], branches[best_branch]->bbox(), bboxes[group]);
        }


        //distribute the rest of braches
        size_t group = nelems[0] < nelems[1] ? group_0 : group_1;
        for (size_t ind = 0; ind < count; ++ind)
        {
            if (!groups[ind])
                groups[ind] = group + 1;
        }
    }

    /*
    PickSeeds pseudo-code:
        PS1 [Calculate ineffiency of grouping entries together] For all pairs of entries
            E1 and E2 a rectangle J is created which includes E1:I and E2:I.
            Calculate
            d = area(J)  area(E1:I)  area(E2:I)
        PS2 [Choose the most wasterful pair] Choose the pair with the largest d
    */
    void pick_seeds(tree_node** branches, size_t count,
                    size_t& group_0, size_t& group_1)
    {
        rect_type       merged;
        coord_type      worst_area = details::lowest<coord_type>::value();

        for (size_t i = 0; i < count; ++i)
        {
            for (size_t j = i + 1; j < count; ++j)
            {
                outersection(branches[i]->bbox(), branches[j]->bbox(), merged);
                coord_type cur_area = area(merged) - area(branches[i]->bbox()) - area(branches[j]->bbox());

                if (cur_area > worst_area)
                {
                    worst_area = cur_area;
                    group_0 = i;
                    group_1 = j;
                }
            }
        }
    }


    /*
    PickNext pseudo-code:
        PN1 [Determine cost of putting each entry in each group] For each entry E which
            is not in a group yet, d1 is calculated. d1 is the area-increase required in
            the covering rectangle of group 1 to include E:I and also d2 for group 2.
        PN2 [Find entry with greatest preference for one 
    */
    size_t pick_next(tree_node** branches, size_t count,
                     size_t const* groups,
                     rect_type const& bbox_0,
                     rect_type const& bbox_1)
    {
        size_t      res = 0;
        coord_type  max_diff = details::lowest<coord_type>::value();

        rect_type   merged;
        for (size_t ind = 0; ind < count; ++ind)
        {
            if (groups[ind] > 0)
            {
                //already picked
                continue;
            }

            coord_type diff = 0;
            outersection(branches[ind]->bbox(), bbox_0, merged);
            diff += area(merged) - area(bbox_0);

            outersection(branches[ind]->bbox(), bbox_1, merged);
            diff -= area(merged) - area(bbox_1);

            if (max_diff < abs(diff))
            {
                max_diff = abs(diff);
                res = ind;
            }
        }

        return res;
    }

    //add child to the node and split node if node is overflowed
    tree_node* add_branch(tree_node* node, tree_node* child)
    {
        assert(node && child && node->level() == child->level() + 1);

        if (node->add_child(child))
        {
            return NULL;
        }

        tree_node* splitted;
        try
        {
            splitted = create_node(node->level());
        }
        catch(...)
        {
            //the only thing we can do here is just to clear child subtree that we can't insert into the parent node
            clean(child);
            throw;
        }

        quadractic_split(node, splitted, child);
        return splitted;
    }


    bool validate_bbox(tree_node const* node)
    {
        if (!node || node->level() == 0)
            return true;


        rect_type rect = recompute_bbox(node);

        if (rect.left() != node->bbox().left() ||
            rect.right() != node->bbox().right() ||
            rect.top() != node->bbox().top() ||
            rect.bottom() != node->bbox().bottom())
            return false;

        return true;
    }


private:
    size_t cached_size_;
    tree_node* root_;

    //tree_node_allocator nalloc_;
    typename Alloc::template rebind<tree_node>::other
        nalloc_;

    R_tree(R_tree const&);
    R_tree& operator=(R_tree const&);
};




//implementation

}//namspace rtree



