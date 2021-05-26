/**
 * @brief Node wrapper
 * 
 * @file node.h
 * @author skyline-nju
 * @date 2018-04-24
 */
#ifndef NODE_H
#define NODE_H
#include "vect.h"
/*************************************************************************//**
 * \brief Unidirectional node wrapper
 * \tparam TPar Template for particles.
 ***************************************************************************/
template <class TPar>
class UniNode : public TPar {
public:
  UniNode() : TPar(), next(nullptr) {}

  UniNode(const TPar &p) : TPar(p), next(nullptr) {}

  UniNode(const double *buf) : TPar(buf), next(nullptr) {}

  template <typename TRan, typename TVec>
  UniNode(TRan &myran, const TVec &l, const TVec &origin)
    : TPar(myran, l, origin), next(nullptr) {}

  void append_at_front(UniNode<TPar> ** head) {
    next = *head;
    *head = this;
  }

  void break_away(UniNode<TPar> *pre_node) const { pre_node->next = next; }

  void break_away(UniNode<TPar> **head, UniNode<TPar> *pre_node) const {
    if (pre_node) {
      pre_node->next = next;
    } else {
      *head = next;
    }
  }
  
  UniNode *next;
};

/**************************************************************************//**
 * \brief Bidirectional node wrapper
 * \tparam TPar Template for particles.
 *****************************************************************************/
template <class TPar>
class BiNode : public TPar {
public:
  BiNode() : TPar(), prev(nullptr), next(nullptr) {}

  BiNode(const TPar &p) : TPar(p), prev(nullptr), next(nullptr) {}

  BiNode(const double *buf) : TPar(buf), prev(nullptr), next(nullptr) {}


  template <typename TRan, typename TVec>
  BiNode(TRan &myran, const TVec &l, const TVec &origin)
    : TPar(myran, l, origin), prev(nullptr), next(nullptr) {}

  void append_at_front(BiNode<TPar> ** head);

  void break_away() const;

  void break_away(BiNode<TPar> **head) const;

  void break_away(BiNode<TPar> **head, BiNode<TPar> *pre_node) const {
    break_away(head);
  }

  BiNode *prev;
  BiNode *next;
};

template <class TPar>
void BiNode<TPar>::append_at_front(BiNode<TPar> ** head) {
  prev = nullptr;
  next = *head;
  if (next) {
    next->prev = this;
  }
  *head = this;
}

template<class TPar>
void BiNode<TPar>::break_away() const {
  prev->next = next;
  if (next) {
    next->prev = prev;
  }
}

template<class TPar>
void BiNode<TPar>::break_away(BiNode<TPar>** head) const {
  if (prev) {
    break_away();
  } else {
    *head = next;
    if (next) {
      next->prev = nullptr;
    }
  }
}

/*************************************************************************//**
 * \brief Visit intra node pairs.
 * \tparam TNode  Template for nodes
 * \tparam BiFunc Template for binary functions: auto(TNode*, TNode*)
 * \param head    Head node of the list, mustn't be a nullptr.
 * \param f_ij    Binary function acting on node i, j
 ****************************************************************************/
template <class TNode, class BiFunc>
void for_each_node_pair(TNode* head, BiFunc f_ij) {
  TNode *node1 = head;
  while (node1->next) {
    TNode *node2 = node1->next;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

/*************************************************************************//**
 * \brief Visit inter node pairs.
 * \tparam TNode    Template for nodes.
 * \tparam BiFunc   Template for binary function: auto(TNode*, TNode*)
 * \param head1     Head node of the first list, mustn't be a nullptr.
 * \param head2     Head node of the second list, mustn't be a nullptr.
 * \param f_ij      Binary function acting on node i, j
 ****************************************************************************/
template <class TNode, class BiFunc>
void for_each_node_pair(TNode* head1, TNode* head2, BiFunc f_ij) {
  TNode *node1 = head1;
  do {
    TNode *node2 = head2;
    do {
      f_ij(node1, node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}

template <class TNode, class TriFunc, class TVec>
void for_each_node_pair(TNode *head1, TNode *head2, const TVec &offset, TriFunc f_ij) {
  TNode *node1 = head1;
  do {
    TNode *node2 = head2;
    do {
      f_ij(node1, node2, offset);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  } while (node1);
}

template<typename TNode>
int count_node(TNode *head) {
  TNode * cur_node = head;
  int count = 0;
  while(cur_node) {
    count++;
    cur_node = cur_node->next;
  }
  return count;
}
#endif