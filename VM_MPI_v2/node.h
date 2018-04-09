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

  void break_away(UniNode<_TPar> *pre_node) const { pre_node->next = next; }
  void break_away(UniNode<_TPar> **head, UniNode<_TPar> *pre_node) const;
  
  UniNode *next;
};

template<class _TPar>
inline void UniNode<_TPar>::break_away(UniNode<_TPar>** head,
                                       UniNode<_TPar>* pre_node) const {
  if (pre_node) {
    pre_node->next = next;
  } else {
    *head = next;
  }
}

template <class _TPar>
class BiNode : public _TPar {
public:
  BiNode() : _TPar(), prev(nullptr), next(nullptr) {}
  BiNode(double x0, double y0, double vx0, double vy0) :
    _TPar(x0, y0, vx0, vy0), prev(nullptr), next(nullptr) {}
  BiNode(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0) :
    _TPar(myran, Lx, Ly, x0, y0), prev(nullptr), next(nullptr) {}
  void append_at_front(BiNode<_TPar> ** head);

  void break_away() const;

  void break_away(BiNode<_TPar> **head) const;

  void break_away(BiNode<_TPar> **head, BiNode<_TPar> *pre_node) const {
    break_away(head);
  }

  BiNode *prev;
  BiNode *next;
};

template <class _TPar>
inline void BiNode<_TPar>::append_at_front(BiNode<_TPar> ** head) {
  prev = nullptr;
  next = *head;
  if (next) {
    next->prev = this;
  }
  *head = this;
}

template<class _TPar>
inline void BiNode<_TPar>::break_away() const {
  prev->next = next;
  if (next) {
    next->prev = prev;
  }
}

template<class _TPar>
inline void BiNode<_TPar>::break_away(BiNode<_TPar>** head) const {
  if (prev) {
    break_away();
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


#endif