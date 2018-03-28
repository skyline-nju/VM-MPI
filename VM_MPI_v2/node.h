#ifndef NODE_H
#define NODE_H
#include "rand.h"

template <class _TPar>
class UniNode : public _TPar {
public:
  UniNode() : _TPar(), next(nullptr) {}
  UniNode(double x0, double y0, double vx0, double vy0) :
    _TPar(x0, y0, vx0, vy0), next(nullptr) {}
  UniNode(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0):
    _TPar(myran, Lx, Ly, x0, y0), next(nullptr) {}

  UniNode *next;
};

template <class _TPar>
class BNode : public _TPar {
public:
  BNode() : Par1(), prev(nullptr), next(nullptr) {}
  BNode(double x0, double y0, double vx0, double vy0) :
    _TPar(x0, y0, vx0, vy0), prev(nullptr), next(nullptr) {}
  BNode(Ran &myran, double Lx, double Ly, double x0 = 0, double y0 = 0) :
    _TPar(myran, Lx, Ly, x0, y0), prev(nullptr), next(nullptr) {}

  BNode* append_front(BNode<_TPar> *head);

  void append_front(BNode<_TPar> **head);

  void break_away(BNode<_TPar> **head) const;

  BNode *prev;
  BNode *next;
};

template<class _TPar>
inline BNode<_TPar> * BNode<_TPar>::append_front(BNode<_TPar> * head) {
  prev = nullptr;
  next = head;
  if (head) {
    head->prev = this;
  }
  return this;
}
template<class _TPar>
inline void BNode<_TPar>::append_front(BNode<_TPar>** head) {
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

#endif
