#ifndef LINKED_LIST_POOL_H
#define LINKED_LIST_POOL_H

#include <vector>
#include <cstdint>
#include <iterator>

/**
 * @brief Pool of singly linked lists.
 *
 * This class is not intended for general-purpose use. It is designed for algorithms that require
 * the ability to split and concatenate linked lists. All data is owned by LinkedListPool.
 *
 * Unlike std::list and std::forward_list, it does not allocate each node separately. Instead,
 * LinkedListPool can reserve memory for multiple lists during construction, using std::vector
 * as the backing container.
 */
template <class T>
class LinkedListPool {
    using IndexType = size_t;

    struct Item {
        IndexType next;
        T value;
    };

public:
    /**
     * @brief List iterator.
     *
     * Iterators remain valid when items are added to the list, although items may be relocated.
     */
    class ListIterator {
        IndexType index = 0;
        LinkedListPool<T>* pool = nullptr;

        ListIterator(IndexType index, LinkedListPool<T>* pool) : index(index), pool(pool) {}

        friend class LinkedListPool<T>;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = size_t;
        using pointer = T*;
        using reference = T&;

        ListIterator() = default;

        reference operator*() { return pool->data[index].value; }
        pointer operator->() { return &pool->data[index].value; }

        ListIterator& operator++() {
            index = pool->data[index].next;
            return *this;
        }

        ListIterator operator++(int) {
            ListIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        bool operator!=(const ListIterator& other) const {
            return index != other.index || pool != other.pool;
        }

        /**
         * @brief Test if iterator points to a valid value.
         */
        explicit operator bool() const { return index != 0; }
    };

    /**
     * @brief Single list within LinkedListPool.
     *
     * A List refers only to a chain of elements. Copying it does not copy any elements. Item data
     * is owned by LinkedListPool.
     *
     * Use LinkedListPool::makeList to create a non-empty list.
     */
    class List {
        IndexType head = 0;
        IndexType tail = 0;

        friend class LinkedListPool;

        List(IndexType head, IndexType tail) : head(head), tail(tail) {}

    public:
        /**
         * @brief Create an empty list.
         */
        List() = default;

        bool isEmpty() const { return head == 0; }
    };

    /**
     * @brief Create a linked list pool with capacity for \a initialCapacity list items.
     * @param initialCapacity Number of elements to preallocate.
     */
    LinkedListPool(size_t initialCapacity) : data(1) {
        data.reserve(initialCapacity + 1); // Reserve space for [0] element
    }

    /**
     * @brief Create a list containing a single item.
     *
     * Does not invalidate any iterators, but may cause item relocation when initialCapacity is
     * exceeded.
     * @param value Value of the element to be inserted in the created list.
     * @return List containing the single value \a value.
     */
    List makeList(const T& value) {
        size_t position = data.size();
        data.push_back(Item{0, value});
        return {position, position};
    }

    /**
     * @brief Split list and return the second half.
     *
     * After performing the operation, the list passed as an argument and the returned list point
     * to the same items. Modifying them will affect both lists.
     *
     * @param list The list to be split.
     * @param head Iterator to the first item in the new list. Needs to be within \a list.
     * @return Returns the suffix of \a list.
     */
    List splitTail(const List& list, const ListIterator& head) {
        return List{head.index, list.tail};
    }

    /**
     * @brief Split list and return the first half.
     *
     * @param list The list to be split.
     * @param end Iterator to the first item that should not be included in the returned list. Needs
     * to be within \a list.
     * @return Returns the prefix of \a list.
     */
    List splitHead(const List& list, const ListIterator& end) {
        if (!end) {
            return list;
        }
        if (end.index == list.head) {
            return {};
        }
        auto last = list.head;
        while (data[last].next != end.index) {
            last = data[last].next;
        }
        data[last].next = 0;
        return List{list.head, last};
    }

    /**
     * @brief Create list iterator from list.
     * @param list
     * @return Iterator pointing to the first item in the list.
     */
    ListIterator head(const List& list) { return iteratorFromIndex(list.head); }

    ListIterator end(const List& list) { return std::next(iteratorFromIndex(list.tail)); }

    List append(const List& head, const List& tail) {
        if (head.isEmpty()) {
            return tail;
        }
        if (tail.isEmpty()) {
            return head;
        }
        List result{head.head, tail.tail};
        data[head.tail].next = tail.head;
        return result;
    }

private:
    ListIterator iteratorFromIndex(IndexType index) { return ListIterator{index, this}; }

    std::vector<Item> data;
};

#endif // LINKED_LIST_POOL_H
