#ifndef MINHEAP_H_
#define MINHEAP_H_
#include <iostream>
#include <vector>
using namespace std;
#include <limits.h>

template<class T, class C>
class MaxHeap {
public:
    MaxHeap() : A(), keys(), heap_size(0) { }
    T &maximum() { 
        if (heap_size) return A[0];
        else {
            cerr << "heap underflow" << endl;
            exit(1);
        }
    }
    T extract_max(C &maxvalue) {
        if (heap_size) {
            T max=A[0];
            maxvalue = keys[0];
            A[0]=A[heap_size-1];
            keys[0]=keys[heap_size-1];
            heap_size--;
            max_heapify(0);
            return max;
        }
        else {
            cerr << "heap underflow" << endl;
            exit(1);
        }
    }
    void decrease_key(int index, C &key) {
        if (key > keys[index]) {
        	cout << "old key = " << keys[index] << " new key = " << key << endl;
            cerr << "new key is smaller than current key" << endl;
            exit(1);
        }
        keys[index]=key;
        while(index>0 && keys[parent(index)]>keys[index]) {
            T temp;
            temp=A[index];
            A[index]=A[parent(index)];
            A[parent(index)]=temp;
            C temp2;
            temp2=keys[index];
            keys[index]=keys[parent(index)];
            keys[parent(index)]=temp2;
            index=parent(index);
        }
    }
    
    void increase_key(int index, C &key) {
        if (key < keys[index]) {
            cerr << "new key is larger than current key" << endl;
            exit(1);
        }
        keys[index]=key;
        while(index <= heap_size-1) {
        	int j = left(index);
        	if(j > heap_size-1) return;
        	if(keys[j] > keys[j+1]) j++; // j contains the smaller children
        	if(j > heap_size-1) return;
        	if(keys[index] <= keys[j]) return;
            T temp;
            temp=A[index];
            A[index]=A[j];
            A[j]=temp;
            C temp2;
            temp2=keys[index];
            keys[index]=keys[j];
            keys[j]=temp2;
            index=j;
        }
    }
    
    void insert(T& x, C &key) {
        heap_size++;
        if (heap_size>A.size()) {
            A.push_back(x);
            keys.push_back(key);
        }
        else {
            A[heap_size-1]=x;
            keys[heap_size-1]=INT_MAX;
        }
        decrease_key(heap_size-1,key);
    }
    
    void update(const T &t, C &key){
    	for(int i=0;i<A.size();i++){
    		if(A[i] == t){
    			A[i] = t;
    			if (key > keys[i]){
//  					cout << "inc old key = " << keys[i] << " new key = " << key << endl;
   					increase_key(i, key);
   				}
    			else{
//   					cout << "dec old key = " << keys[i] << " new key = " << key << endl;
    				decrease_key(i, key);
    			}
    			return;	
    		}    			
    	}
    	return;
    }
    bool empty() const{
    	return heap_size == 0;
    }
    
    int size() const { return heap_size; }
private:
    int parent(int i) { return (i-1)/2; }
    int left(int i) { return 2*i+1; }
    int right(int i) { return 2*i+2; }
    void build_max_heap() {
        heap_size=A.size();
        for(int i=heap_size/2;i<=0;--i)
            max_heapify(i);
    }
    vector<T> A;
    vector<C> keys;
    int heap_size;
    void max_heapify(int i) {
        int l = left(i);
        int r = right(i);
        int smallest;
        if (l<heap_size && keys[l] < keys[i])
            smallest = l;
        else
            smallest = i;
        if (r<heap_size && keys[r] < keys[smallest])
            smallest = r;
        if (smallest != i) {
            T temp;
            temp=A[i];
            A[i]=A[smallest];
            A[smallest]=temp;

            C temp2;
            temp2=keys[i];
            keys[i]=keys[smallest];
            keys[smallest]=temp2;
    
            min_heapify(smallest);
        }
    }
};

#endif /*MinHeap_H_*/
