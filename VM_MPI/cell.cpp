#include "cell.h"
using namespace std;

void Cell::interact() {
}

void Cell::interact(Cell *c2) {

}

void Cell::mark_as_invalid(stack <unsigned int> &empty_pos) {
  if (head) {
    Node *curNode = head;
    do {
      curNode->valid = false;
      empty_pos.push(curNode->par_idx);
      curNode->next = curNode;
    } while (curNode);
    head = nullptr;
  }
}
void Cell::move(Ran* myran, double eta, Cell *cell, double Lx, double Ly,
                double Ly_l, double Ly_h) {
  if (head) {
    Node *curNode = head;
    Node *preNode = nullptr;
    do {
      double noise = eta * 2 * PI * (myran->doub() - 0.5); //disorder free
      bool out_range;
      curNode->update_coor(noise, Lx, Ly, Ly_l, Ly_h, out_range);
      if (out_range) {
        if (preNode) {
          preNode->next = curNode->next;
          curNode->next = cell[curNode->cell_idx].head;
          cell[curNode->cell_idx].head = curNode;
          curNode = preNode->next;
        }
        else {
          head->next = curNode->next;
          curNode->next = cell[curNode->cell_idx].head;
          cell[curNode->cell_idx].head = curNode;
          curNode = head->next;
        }
      } else {
        preNode = curNode;
        curNode = curNode->next;
      }
    } while (curNode);
   }
}