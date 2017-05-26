#include "grid.h"

using namespace std;

int Grid::mx;
int Grid::my;
int Grid::mm;
double Grid::l0 = 1.0;

Grid * Grid::ini(double Lx, double Ly) {
  mx = int(Lx / l0);
  my = int(Ly / l0);
  mm = mx * my;
  Grid *cell = new Grid[mm];
  for (int i = 0; i < mm; i++) {
    cell[i].head = nullptr;
  }
  return cell;
}

void Grid::cell_cell() {
  Node *node1 = head;
  Node *node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      node1->align(node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

void Grid::cell_cell(Grid* grid) {
  if (grid->head) {
    Node* node1 = head;
    Node* node2;
    do {
      node2 = grid->head;
      do {
        node1->align(node2);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

void Grid::cell_cell(Grid* grid, double a, double b) {
  if (grid->head) {
    Node* node1 = head;
    Node* node2;
    do {
      node2 = grid->head;
      do {
        node1->align(node2, a, b);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

void Grid::all_pairs(Grid * cell) {
  int i, j;
  Grid* p = cell;
  for (j = 0; j <= my - 2; j++) {
    if (p->head) {
      p->cell_cell();
      p->cell_cell(p + 1);
      p->cell_cell(p + mx + mx - 1, -Node::Lx, 0);
      p->cell_cell(p + mx);
      p->cell_cell(p + mx + 1);
    }
    p++;
    for (i = 1; i <= mx - 2; i++) {
      if (p->head) {
        p->cell_cell();
        p->cell_cell(p + 1);
        p->cell_cell(p + mx - 1);
        p->cell_cell(p + mx);
        p->cell_cell(p + mx + 1);
      }
      p++;
    }
    if (p->head) {
      p->cell_cell();
      p->cell_cell(p - mx + 1, Node::Lx, 0);
      p->cell_cell(p + 1, Node::Lx, 0);
      p->cell_cell(p + mx);
      p->cell_cell(p + mx - 1);
    }
    p++;
  }
  if (p->head) {
    p->cell_cell();
    p->cell_cell(p + 1);
    p->cell_cell(cell, 0, Node::Ly);
    p->cell_cell(cell + 1, 0, Node::Ly);
    p->cell_cell(cell + mx - 1, -Node::Lx, Node::Ly);
  }
  p++;
  for (i = 1; i <= mx - 2; i++) {
    if (p->head) {
      p->cell_cell();
      p->cell_cell(p + 1);
      p->cell_cell(cell + i - 1, 0, Node::Ly);
      p->cell_cell(cell + i, 0, Node::Ly);
      p->cell_cell(cell + i + 1, 0, Node::Ly);
    }
    p++;
  }
  if (p->head) {
    p->cell_cell();
    p->cell_cell(p - mx + 1, Node::Lx, 0);
    p->cell_cell(cell + mx - 2, 0, Node::Ly);
    p->cell_cell(cell + mx - 1, 0, Node::Ly);
    p->cell_cell(cell, Node::Lx, Node::Ly);
  }
}

void Grid::link_nodes(Grid * cell, Node * node) {
  for (int i = 0; i < Node::N; i++) {
    int col = int(node[i].pdata.x);
    if (col >= mx)
      col -= mx;
    else if (col < 0)
      col += mx;
    int row = int(node[i].pdata.y);
    if (row >= my)
      row -= my;
    else if (row < 0)
      row += my;
    int j = node[i].cell_idx = col + mx * row;
    node[i].next = cell[j].head;
    cell[j].head = &node[i];
  }
}

void Grid::refresh(Grid *cell, Node *node) {
  for (int i = 0; i < mm; i++) {
    cell[i].head = nullptr;
  }
  Grid::link_nodes(cell, node);
}
