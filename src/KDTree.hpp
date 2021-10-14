// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
  size_t dimension_;
  size_t size_;
  std::vector<value_type> data;

  struct Node {
        Point<N> point;
        Node *left;
        Node *right;
        int level;  // level of the node in the tree, starts at 0 for the root
        ElemType value;
        Node(const Point<N>& _pt, int _level, const ElemType& _value=ElemType()){
          point = _pt;
          left = nullptr;
          right = nullptr;
          level = _level;
          value = _value;
        }
    };

    Node* root_;

    Node* clone(Node* root);

};

template <std::size_t N, typename ElemType>
typename KDTree<N, ElemType>::Node* KDTree<N, ElemType>::clone(typename KDTree<N, ElemType>::Node* root) {
    if (root == NULL){
      return NULL;
    }
    Node* newRoot = new Node(*root);
    newRoot->left = clone(root->left);
    newRoot->right = clone(root->right);
    return newRoot;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
	dimension_ = N;
	size_ = 0;
  root_ = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // (me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // (me): Fill this in.
  size_ = rhs.size_;
  dimension_ = rhs.dimension_;
  data = rhs.data;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    size_ = rhs.size_;
    dimension_ = rhs.dimension_;
    data = rhs.data;
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  // (me): Fill this in.
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  // (me): Fill this in.
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  // (me): Fill this in.
  if(size_ == 0) {
	 return true;
  }else{
  	return false;
  }
	
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  // (me): Fill this in.
  for (size_t i = 0; i < data.size(); ++i) {
  	if(data[i].first == pt){
  	     return true;
  	}
  }
  return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  // (me): Fill this in.
  value_type temp(pt,value);
  if(contains(pt)){
    data[at(pt)].first = pt;
    data[at(pt)].second = value;
  }else{
    data.push_back(temp);
    size_ = data.size();
  }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  // (me): Fill this in.
  if(contains(pt)){
    insert(pt,data[at(pt)].second);
  }else{
    insert(pt,data.size());
  }
  return at(pt);
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  // (me): Fill this in.
  
  for (ElemType i = 0; i < data.size(); ++i) {
  	if(data[i].first == pt){
  	     return data[i].second;
  	}
  }
  throw std::out_of_range ("No existe");
  
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  // (me): Fill this in.
  for (ElemType i = 0; i < data.size(); ++i) {
    if(data[i].first == pt){
         return data[i].second;
    }
  }
  throw std::out_of_range ("No existe");
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // (me): Fill this in.
  if(empty()){
    return ElemType();
  }
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key, size_t k) const {
  // (me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// (me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
