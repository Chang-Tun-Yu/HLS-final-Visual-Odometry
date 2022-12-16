#include <iostream>
#include <iomanip>
#include <cstdlib>
//#include "ap_int.h"
using namespace std;

#define MaxSize 500
//Vector Class
//template<class data_T, class flag_T>
template<class data_T>
class Vector {
protected:
	data_T storage[MaxSize];
	bool empty();
	bool full();
	int size;

private:
public:
	void     clear();
	void     push_back(data_T element);
	data_T   fetch(int idx);
	data_T   front();
	data_T   back();
	data_T   *begin();
	data_T   *end();

};
template<class data_T>
void   Vector<data_T>::clear() {
	size = 0;
}
template<class data_T>
void Vector<data_T>::push_back(data_T element) {
	bool flag = !full();
	if (flag)
	{
		storage[size] = element;
		size = size + 1;
	}
}
template<class data_T>
data_T Vector<data_T>::front() {
	bool flag = !empty();
	data_T data;
	if (flag)
	{
		data = storage[0];
	}
	return data;
}
template<class data_T>
data_T Vector<data_T>::fetch(int idx) {
	bool flag = !((empty()) && (idx < size - 1));
	data_T data;
	if (flag)
	{
		data = storage[idx];
	}
	return data;
}
template<class data_T>
data_T Vector<data_T>::back() {
	bool flag = !empty();
	data_T data;
	if (flag)
	{
		data = storage[size - 1];
	}
	return data;
}
template<class data_T>
data_T* Vector<data_T>::begin() {
	bool flag = !empty();
	data_T* iter = 0;
	if (flag)
	{
		iter = storage;
	}
	return iter;
}
template<class data_T>
data_T* Vector<data_T>::end() {
	bool flag = !empty();
	data_T* iter = 0;
	if (flag)
	{
		iter = (storage + size);
	}
	return iter;
}
template<class data_T>
bool Vector<data_T>::empty() {
	return (size == 0);
}
template<class data_T>
bool Vector<data_T>::full() {
	return (size == MaxSize);
}
