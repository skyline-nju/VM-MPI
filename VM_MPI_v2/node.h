#ifndef NODE_H
#define NODE_H
#include "rand.h"
#include <list>
#include <vector>

template <class _TPar>
class UniNode : public _TPar {
public:
  UniNode() : _TPar(), next(nullptr) {}
  UniNode(double x0, double y0, double vx0, double vy0) :
    _TPar(x0, y0, vx0, vy0), next(nullptr) {}
  UniNode(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0):
    _TPar(myran, Lx, Ly, x0, y0), next(nullptr) {}

  void append_at_front(UniNode<_TPar> ** head) {
    next = *head;
    *head = this;
  }

  void break_away(UniNode<_TPar> **head, UniNode<_TPar> *pre_node);
  
  UniNode *next;
};

template<class _TPar>
inline void UniNode<_TPar>::break_away(UniNode<_TPar>** head,
                                       UniNode<_TPar>* pre_node) {
  if (pre_node) {
    pre_node->next = next;
  } else {
    *head = next;
  }
}

template <class _TPar>
class BNode : public _TPar {
public:
  BNode() : Par1(), prev(nullptr), next(nullptr) {}
  BNode(double x0, double y0, double vx0, double vy0) :
    _TPar(x0, y0, vx0, vy0), prev(nullptr), next(nullptr) {}
  BNode(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0) :
    _TPar(myran, Lx, Ly, x0, y0), prev(nullptr), next(nullptr) {}
  void append_at_front(BNode<_TPar> ** head);

  void break_away(BNode<_TPar> **head) const;

  void break_away(BNode<_TPar> **head, BNode<_TPar> *pre_node) {
    break_away(head);
  }

  BNode *prev;
  BNode *next;
};

template <class _TPar>
inline void BNode<_TPar>::append_at_front(BNode<_TPar> ** head) {
  prev = nullptr;
  next = *head;
  if (next) {
    next->prev = this;
  }
  *head = this;
}


template<class _TPar>
inline void BNode<_TPar>::break_away(BNode<_TPar>** head) const {
  if (prev) {
    prev->next = next;
    if (next) {
      next->prev = prev;
    }
  } else {
    *head = next;
    if (next) {
      next->prev = nullptr;
    }
  }
}

template <class _TNode, class BiFunc>
void for_each_node_pair(_TNode* head, BiFunc f_ij) {
  _TNode *node1 = head;
  while (node1->next) {
    _TNode *node2 = node1->next;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template <class _TPar, class BiFunc>
void for_each_node_pair(const std::list<_TPar *> &cl, BiFunc f_ij) {
  auto end2 = cl.cend();
  auto end1 = std::prev(end2);
  for (auto it1 = cl.cbegin(); it1 != end1; ++it1) {
    for (auto it2 = std::next(it1); it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template <class _TNode, class BiFunc>
void for_each_node_pair(_TNode* head1, _TNode* head2, BiFunc f_ij) {
  _TNode *node1 = head1;
  do {
    _TNode *node2 = head2;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}

template <class _TPar, class BiFunc>
void for_each_node_pair(const std::list<_TPar *> &cl1,
                        const std::list<_TPar *> &cl2, BiFunc f_ij) {
  auto end1 = cl1.cend();
  auto end2 = cl2.cend();
  auto beg2 = cl2.cbegin();
  for (auto it1 = cl1.cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template <class _TNode>
void del_list(std::vector<_TNode *> &list_arr) {
  auto end = list_arr.end();
  for (auto it = list_arr.begin(); it != end; ++it) {
    *it = nullptr;
  }
}

template <class _TPar>
void del_list(std::vector<std::list<_TPar *>> &list_arr) {
  auto end = list_arr.end();
  for (auto it = list_arr.begin(); it != end; ++it) {
    (*it).clear();
  }
}

template <class T>
inline bool not_empty(const T * ptr) {
  return ptr != nullptr;
}

template <class T>
inline bool not_empty(const std::list<T *> &l) {
  return !l.empty();
}

template <class T>
inline void merge_list(UniNode<T> **head, UniNode<T> **tail, UniNode<T> **head_tmp) {
  if (*head_tmp) {
    if (*head) {
      (*tail)->next = *head_tmp;
    } else {
      *head = *head_tmp;
    }
    *head_tmp = nullptr;
  }
  *tail = nullptr;
}

template <class T>
inline void merge_list(BNode<T> **head, BNode<T> **tail, BNode<T> **head_tmp) {
  if (*head_tmp) {
    if (*head) {
      (*tail)->next = *head_tmp;
      (*head_tmp)->prev = *tail; 
    } else {
      *head = *head_tmp;
    }
    *head_tmp = nullptr;
  }
  *tail = nullptr;
}

#endif


