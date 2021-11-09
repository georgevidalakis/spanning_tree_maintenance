#ifndef MST_SUBSET_ITERATOR_H_
#define MST_SUBSET_ITERATOR_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#include <set>

template<typename T> class SubsetIterator {
 public:
  SubsetIterator() = delete;
  
  SubsetIterator(std::vector<T> set_items) : set_items_(set_items), set_size_(set_items_.size()) {
    belongs_to_subset_ = (bool*)malloc(set_size_ * sizeof(bool));
    if (set_size_ && belongs_to_subset_ == NULL) {
      std::cerr << "[ERROR]: SubsetIterator could not allocate the required memory.\n";
      exit(1);
    }
    for (int item_idx = 0; item_idx < set_size_; ++item_idx) {
      belongs_to_subset_[item_idx] = false;
    }
    reached_end_of_iterator_ = !set_size_;
  }
  
  ~SubsetIterator() {
    free(belongs_to_subset_);
  }
  
  bool HasReachedEndOfIterator() const {
    return reached_end_of_iterator_;
  }
  
  std::set<T> GetSubset() {
    if (reached_end_of_iterator_) {
      std::cerr << "[ERROR]: SubsetIterator's GetSubset method was called while the end of the iterator had been reached.\n";
      exit(1);
    }
    std::set<T> subset;
    for (int item_idx = 0; item_idx < set_size_; ++item_idx) {
      if (belongs_to_subset_[item_idx]) {
        subset.insert(set_items_[item_idx]);
      }
    }
    NextSubset();
    return subset;
  }

 private:
  void NextSubset() {
    if (reached_end_of_iterator_) {
      std::cerr << "[ERROR]: SubsetIterator's NextSubset method was called while the end of the iterator had been reached.\n";
      exit(1);
    }
    reached_end_of_iterator_ = true;
    for (int item_idx = 0; item_idx < set_size_; ++item_idx) {
      if (belongs_to_subset_[item_idx]) {
        belongs_to_subset_[item_idx] = false;
      }
      else {
        belongs_to_subset_[item_idx] = true;
        reached_end_of_iterator_ = false;
        break;
      }
    }
  }
  
  std::vector<T> set_items_;
  const int set_size_;
  bool *belongs_to_subset_;
  bool reached_end_of_iterator_;
};

#endif  // MST_SUBSET_ITERATOR_H_
