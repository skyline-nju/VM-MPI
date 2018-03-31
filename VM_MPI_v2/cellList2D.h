#ifndef CELLLIST2D_H
#define CELLLIST2D_H
#include <list>
#include <vector>
#include <iostream>
#include <typeinfo>
#include "vect.h"
#include "comn.h"
#include "node.h"

template <typename _TCell>
class Cell_list_2 {
public:
	Cell_list_2(double Lx, double Ly, double x0 = 0, double y0 = 0, double r_cut = 1,
		bool comm_x = false, bool comm_y = false);

	template <typename T>
	int get_ic(const T &p_arr) const;

	template <typename T>
	void create(std::vector<T> &p_arr) { create_lists(p_arr, cell); }

	template <typename T>
	void recreate(std::vector<T> &p_arr) { del_list(cell); create(p_arr); }

	template <typename UniFunc>
	void update(UniFunc move) { update_lists(cell, move); }

	template <typename BiFunc>
	void for_each_pair(BiFunc f_ij, int row_beg = 0, int row_end = -1, int col_beg = 0, int col_end = -1);

private:
	template <typename _TPar>
	void create_lists(std::vector<_TPar> &p_arr, std::vector<std::list<_TPar*>> &cl);

	template <typename _TNode>
	void create_lists(std::vector<_TNode> &p_arr, std::vector<_TNode*> &cl);

	template <typename _TPar, typename UniFunc>
	void update_lists(std::vector<std::list<_TPar*>> &cl, UniFunc move);

	template <typename _TNode, typename UniFunc>
	void update_lists(std::vector<_TNode*> &cl, UniFunc move);


protected:
	Vec_2<double> origin;
	Vec_2<double> l_box;
	Vec_2<double> inverse_lc;
	Vec_2<double> l_cell;
	Vec_2<int> n_bins;
	int ncells;
	std::vector<_TCell> cell;
};

template <typename _TCell>
Cell_list_2<_TCell>::Cell_list_2(double Lx, double Ly, double x0, double y0,
	double r_cut, bool comm_x, bool comm_y) {
	if (comm_x == false && comm_y == false) {
		origin.x = 0;
		origin.y = 0;
		n_bins.x = int(Lx / r_cut);
		n_bins.y = int(Ly / r_cut);
		l_cell.x = Lx / n_bins.x;
		l_cell.y = Ly / n_bins.y;
		l_box.x = Lx;
		l_box.y = Ly;
	} else if (comm_x == false && comm_y == true) {
		n_bins.x = int(Lx / r_cut);
		n_bins.y = int(Ly / r_cut) + 2;
		l_cell.x = Lx / n_bins.x;
		l_cell.y = Ly / (n_bins.y - 2);
		l_box.x = Lx;
		l_box.y = Ly + 2 * l_cell.y;
		origin.x = x0;
		origin.y = y0 - l_cell.y;
	} else if (comm_x == true && comm_y == true) {
		n_bins.x = int(Lx / r_cut) + 2;
		n_bins.y = int(Ly / r_cut) + 2;
		l_cell.x = Lx / (n_bins.x - 2);
		l_cell.y = Ly / (n_bins.y - 2);
		l_box.x = Lx + 2 * l_cell.x;
		l_box.y = Ly + 2 * l_cell.y;
		origin.x = x0 - l_cell.x;
		origin.y = y0 - l_cell.y;
	}
	inverse_lc.x = 1 / l_cell.x;
	inverse_lc.y = 1 / l_cell.y;
	ncells = n_bins.x * n_bins.y;
	cell.reserve(ncells);
	for (int i = 0; i < ncells; i++) {
		cell.emplace_back();
	}

	std::cout << "L_Box = " << l_box.x << "\t" << l_box.y << "\n"
		<< "L_cell = " << l_cell.x << "\t" << l_cell.y << "\nx0 = "
		<< origin.x << "\t y0 = " << origin.y << "\nnx = " << n_bins.x
		<< "\tny = " << n_bins.y << "\n";
}

template<typename _TCell>
template<typename T>
inline int Cell_list_2<_TCell>::get_ic(const T & p) const {
	return int((p.x - origin.x) * inverse_lc.x)
		+ int((p.y - origin.y) * inverse_lc.y) * n_bins.x;
}

template<typename _TCell>
template<typename _TPar>
void Cell_list_2<_TCell>::create_lists(std::vector<_TPar>& p_arr,
	std::vector<std::list<_TPar*>>& cl) {
	auto end = p_arr.end();
	for (auto it = p_arr.begin(); it != end; ++it) {
		int ic = get_ic(*it);
		cl[ic].push_back(&(*it));
	}
}

template<typename _TCell>
template <typename _TNode>
void Cell_list_2<_TCell>::create_lists(std::vector<_TNode> &p_arr, std::vector<_TNode*> &cl) {
	auto end = p_arr.end();
	for (auto it = p_arr.begin(); it != end; ++it) {
		int ic = get_ic(*it);
		(*it).append_at_front(&cl[ic]);
	}
}

template<typename _TCell>
template<typename _TPar, typename UniFunc>
void Cell_list_2<_TCell>::update_lists(std::vector<std::list<_TPar*>>& cl, UniFunc move) {
	int *list_len = new int[ncells];
	for (int ic = 0; ic < ncells; ic++) {
		list_len[ic] = cl[ic].size();
	}
	for (int ic = 0; ic < ncells; ic++) {
		int depth = list_len[ic];
		if (depth) {
			int it_count = 0;
			for (auto it = cl[ic].begin(); it_count < depth; it_count++) {
				_TPar *ptr = *it;
				move(ptr);
				int ic_new = get_ic(*ptr);
				if (ic_new == ic) {
					++it;
				} else {
					it = cl[ic].erase(it);
					cl[ic_new].push_back(ptr);
				}
			}
		}
	}
	delete[] list_len;
}

template<typename _TCell>
template<typename _TNode, typename UniFunc>
void Cell_list_2<_TCell>::update_lists(std::vector<_TNode*>& cl, UniFunc move) {
	std::vector<_TNode*> tail(ncells, nullptr);
	std::vector<_TNode*> head_tmp(ncells, nullptr);
	for (int ic = 0; ic < ncells; ic++) {
		if (cl[ic]) {
			auto cur_node = cl[ic];
			_TNode * pre_node = nullptr;
			do {
				move(cur_node);
				int ic_new = get_ic(*cur_node);
				if (ic_new != ic) {
					cur_node->break_away(&cl[ic], pre_node);
					auto tmp_node = cur_node;
					cur_node = cur_node->next;
					tmp_node->append_at_front(&head_tmp[ic_new]);
				} else {
					pre_node = cur_node;
					cur_node = cur_node->next;
				}
			} while (cur_node);
			tail[ic] = pre_node;
		}
	}
	for (int ic = 0; ic < ncells; ic++) {
		merge_list(&cl[ic], &tail[ic], &head_tmp[ic]);
	}
}

template<typename _TCell>
template<typename BiFunc>
void Cell_list_2<_TCell>::for_each_pair(BiFunc f_ij, int row_beg, int row_end, int col_beg, int col_end) {
	if (row_end == -1)
		row_end = n_bins.y;
	if (col_end == -1)
		col_end = n_bins.x;
	for (int row = row_beg; row < row_end; row++) {
		int row_upper = row + 1;
		if (row_upper == n_bins.y)
			row_upper = 0;
		int row_times_ncols = row * n_bins.x;
		int row_upper_times_ncols = row_upper * n_bins.x;
		for (int col = col_beg; col < col_end; col++) {
			int ic0 = col + row_times_ncols;
			int col_right = col + 1;
			if (col_right == n_bins.x)
				col_right = 0;
			int ic1 = col_right + row_times_ncols;
			int ic2 = col + row_upper_times_ncols;
			if (not_empty(cell[ic0])) {
				for_each_node_pair(cell[ic0], f_ij);
				int ic3 = col_right + row_upper_times_ncols;
				bool flag_c1_c2 = true;
				if (not_empty(cell[ic1]))
					for_each_node_pair(cell[ic0], cell[ic1], f_ij);
				else
					flag_c1_c2 = false;
				if (not_empty(cell[ic2]))
					for_each_node_pair(cell[ic0], cell[ic2], f_ij);
				else
					flag_c1_c2 = false;
				if (flag_c1_c2)
					for_each_node_pair(cell[ic1], cell[ic2], f_ij);
				if (not_empty(cell[ic3]))
					for_each_node_pair(cell[ic0], cell[ic3], f_ij);
			} else if (not_empty(cell[ic1]) && not_empty(cell[ic2])) {
				for_each_node_pair(cell[ic1], cell[ic2], f_ij);
			}
		}
	}
}

#endif
